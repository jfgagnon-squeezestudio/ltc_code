#include "brdf_hair.h"

#define M_PI        3.14159265358979323846f
#define M_2PI       6.28318530717958647692f
#define M_1_2PI     0.159154943091895335768f

#define M_SQRT_PI_8	0.626657069f
#define M_LN_2PI	1.8378770664093454f // ln(2*pi)
#define P_MAX		3

void Assert(bool test)
{
	if (!test)
	{
		// Just for quick debug/test
		int* i = nullptr;
		*i = 0;
	}
}


struct PrincipledHairBsdf
{
	float h;
	float eta;
	float beta_m;
	float beta_n;
	float alpha;
	vec3 sigma_a;
	vec3 tangent;

	float gammaO;
	float v[P_MAX + 1];
	float s;

	// TODO: THIS SHOULD GO CPU SIDE
	float sin2kAlpha[3];
	float cos2kAlpha[3];
};

// Declare some usefull math functions. Some are to tolerate out of range values due to floating
// point round-off error, some are convenience

float SafeASin(float angle)
{
	angle = clamp(angle, -1.0f, 1.0f);
	return sin(angle);
}

float SafeSqrt(float value)
{
	value = max(0.0f, value);
	return sqrt(value);
}

float Sqr(float value)
{
	return value * value;
}


// Hair bsdf proper

void PrincipledHairBsdfDecodeGBuffer(float h, vec4 baseColorBuffer, vec4 normalBuffer, vec4 materialBuffer, PrincipledHairBsdf& bsdf)
{
	bsdf.h = h; // TODO: find a way to compute this from gbuffer
	bsdf.eta = normalBuffer.a;
	bsdf.beta_m = materialBuffer.r;
	bsdf.beta_n = materialBuffer.g;
	bsdf.alpha = materialBuffer.a;
	bsdf.sigma_a = vec3(baseColorBuffer);
	bsdf.tangent = vec3(normalBuffer);

	bsdf.gammaO = SafeASin(bsdf.h);

	// Compute longitudinal variance
	bsdf.v[0] = Sqr(0.726 * bsdf.beta_m + 0.812 * Sqr(bsdf.beta_m) + 3.7 * pow(bsdf.beta_m, 20));
	bsdf.v[1] = 0.25 * bsdf.v[0];
	bsdf.v[2] = 4.0 * bsdf.v[0];
	for (int p = 3; p <= P_MAX; ++p)
		bsdf.v[p] = bsdf.v[2];

	// Compute azimuthal logistic scale factor from beta_n
	bsdf.s = M_SQRT_PI_8 * (0.265f * bsdf.beta_n + 1.194f * Sqr(bsdf.beta_n) + 5.372f * pow(bsdf.beta_n, 22));

	// Compute alpha terms for hair scales (TODO: THIS SHOULD GO CPU SIDE)
	bsdf.sin2kAlpha[0] = sin(bsdf.alpha);
	bsdf.cos2kAlpha[0] = SafeSqrt(1 - Sqr(bsdf.sin2kAlpha[0]));
	for (int i = 1; i < 3; ++i)
	{
		bsdf.sin2kAlpha[i] = 2.0 * bsdf.cos2kAlpha[i - 1] * bsdf.sin2kAlpha[i - 1];
		bsdf.cos2kAlpha[i] = Sqr(bsdf.cos2kAlpha[i - 1]) - Sqr(bsdf.sin2kAlpha[i - 1]);
	}
}


// Numerical approximation to the Bessel function of the first kind.
float I0(float x)
{
	float result = 1;
	float factorialPow2 = 1.0;
	float fourPowI = 4.0;
	float xPow2 = Sqr(x);
	float xPow2I = xPow2;

	// See bsdf_hair_principled.h for an alternative implementation
	for (int i = 1; i <= 10; ++i)
	{
		factorialPow2 *= i * i;
		result += xPow2I / (fourPowI * factorialPow2);
		xPow2I *= xPow2;
		fourPowI *= 4.0;
	}

	return result;
}

// Logarithm of the Bessel function of the first kind.
float LogI0(float x)
{
	if (x > 12.0f)
	{
		// log(1/x) == -log(x) iff x > 0.
		return x + 0.5f * (1.f / (8.0f * x) - M_LN_2PI - log(x));
	}
	else
	{
		return log(I0(x));
	}
}

// Change azimutal angle
float Phi(int p, float gammaO, float gammaT)
{
	return 2 * p * gammaT - 2 * gammaO + p * M_PI;
}

// Logistic dsitribution. Used to compute glossy scattering for azimutal events
float Logistic(float x, float s)
{
	// Used abs() for numerical stability. See pbrt.
	x = abs(x);
	return exp(-x / s) / (s * Sqr(1 + exp(-x / s)));
}

float LogisticCDF(float x, float s)
{
	return 1 / (1 + exp(-x / s));
}

float TrimmedLogistic(float x, float s, float a, float b)
{
	return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}

// BxDF Utility Functions
float FrDielectric(float cosThetaI, float etaI, float etaT)
{
	cosThetaI = clamp(cosThetaI, -1.0f, 1.0f);

	// We do not handle wich side ray is coming
	// etaI and etaT mut be set properly by caller
	cosThetaI = abs(cosThetaI);

	// TODO: verify it is acceptable to remove this test
	////// Potentially swap indices of refraction
	////bool entering = cosThetaI > 0.f;
	////if (!entering)
	////{
	////	float tmp = etaI;
	////	etaI = etaT;
	////	etaT = tmp;
	////	cosThetaI = abs(cosThetaI);
	////}

	// Compute cosThetaT using Snell's law
	float sinThetaI = sqrt(max(0.0, 1.0 - cosThetaI * cosThetaI));
	float sinThetaT = etaI / etaT * sinThetaI;

	// Handle total internal reflection
	if (sinThetaT >= 1) return 1.0;
	float cosThetaT = sqrt(max(0.0, 1.0 - sinThetaT * sinThetaT));
	float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) / ((etaT * cosThetaI) + (etaI * cosThetaT));
	float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) / ((etaI * cosThetaI) + (etaT * cosThetaT));
	return (Rparl * Rparl + Rperp * Rperp) * 0.5;
}

// Longitudinal scattering
float Mp(float cosThetaI, float cosThetaO, float sinThetaI, float sinThetaO, float v)
{
	float a = cosThetaI * cosThetaO / v;
	float b = sinThetaI * sinThetaO / v;
	float mp = (v <= .1) ?
		(exp(LogI0(a) - b - 1 / v + 0.6931f + log(1 / (2 * v)))) :
		(exp(-b) * I0(a)) / (sinh(1 / v) * 2 * v);
	return mp;
}

// Absorption
void Ap(float cosThetaO, float eta, float h, vec3 T, vec3* ap/*[P_MAX + 1]*/)
{
	// Compute p = 0 attenuation at initial cylinder intersection
	float cosGammaO = SafeSqrt(1 - h * h);
	float cosTheta = cosThetaO * cosGammaO;

	float f = FrDielectric(cosTheta, 1.0, eta);
	ap[0] = vec3(f);

	// Compute p = 1 attenuation term
	ap[1] = Sqr(1 - f) * T;

	// Compute attenuation terms up to p = P_MAX
	for (int p = 2; p < P_MAX; ++p)
		ap[p] = ap[p - 1] * T * f;

	// Compute attenuation term accounting for remaining orders of scattering
	ap[P_MAX] = ap[P_MAX - 1] * f * T / (vec3(1.f) - T * f);
}

// Azimutal scattering
float Np(float phi, int p, float s, float gammaO, float gammaT)
{
	float dphi = phi - Phi(p, gammaO, gammaT);
	// Remap dphi to[?M_PI, M_PI]
	while (dphi > M_PI) dphi -= M_2PI;
	while (dphi < -M_PI) dphi += M_2PI;
	return TrimmedLogistic(dphi, s, -M_PI, M_PI);
}

vec3 EvalPrincipledHairBsdf(vec3 wo, vec3 wi, PrincipledHairBsdf bsdf)
{
	// bsdf = sum( Mp(theta_o, theta_i) * Ap(omega_o) * Np(phi) )
	//        ---------------------------------------------------
	//                       | cos(theta_i) |
	// for p = [0, 2]
	// Where:
	//     Mp = longitudinal scattering
	//     Ap = absorption
	//     Np = azimutal scattering
	//     p  = number of times internal refraction happens (3 our case)
	// Refer to https://tinyurl.com/yypps93u for a detailed diagram

	// Compute hair coordinate system terms related to wo|i
	float sinThetaO = wo.x;
	float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
	float phiO = atan(wo.z, wo.y);

	float sinThetaI = wi.x;
	float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
	float phiI = atan(wi.z, wi.y);

	// Compute cos(theta_t) for refracted ray
	float sinThetaT = sinThetaO / bsdf.eta;
	float cosThetaT = SafeSqrt(1.0 - Sqr(sinThetaT));

	// Compute gamma_t for refracted ray
	float etap = sqrt(bsdf.eta * bsdf.eta - Sqr(sinThetaO)) / cosThetaO;
	float sinGammaT = bsdf.h / etap;
	float cosGammaT = SafeSqrt(1.0 - Sqr(sinGammaT));
	float gammaT = SafeASin(sinGammaT);

	// Compute the transmittance T of a single path through the cylinder
	vec3 T = exp(-bsdf.sigma_a * (2.0f * cosGammaT / cosThetaT));

	// Evaluate hair BSDF
	vec3 ap[P_MAX + 1];
	Ap(cosThetaO, bsdf.eta, bsdf.h, T, ap);

	float phi = phiI - phiO;
	vec3 fsum = vec3(0.0);
	for (int p = 0; p < P_MAX; ++p)
	{
		// Compute sin(theta_i) and cos(theta_i) terms accounting for scales offset
		float sinThetaIp, cosThetaIp;
		if (p == 0)
		{
			// R -> theta_i += 2*alpha
			sinThetaIp = sinThetaI * bsdf.cos2kAlpha[1] + cosThetaI * bsdf.sin2kAlpha[1];
			cosThetaIp = cosThetaI * bsdf.cos2kAlpha[1] - sinThetaI * bsdf.sin2kAlpha[1];
		}
		else if (p == 1)
		{
			// TT -> theta_i -= alpha
			sinThetaIp = sinThetaI * bsdf.cos2kAlpha[0] - cosThetaI * bsdf.sin2kAlpha[0];
			cosThetaIp = cosThetaI * bsdf.cos2kAlpha[0] + sinThetaI * bsdf.sin2kAlpha[0];
		}
		else if (p == 2)
		{
			// TRT -> theta_i -= 4*alpha
			sinThetaIp = sinThetaI * bsdf.cos2kAlpha[2] - cosThetaI * bsdf.sin2kAlpha[2];
			cosThetaIp = cosThetaI * bsdf.cos2kAlpha[2] + sinThetaI * bsdf.sin2kAlpha[2];
		}

		// Handle out-of-range cos(theta_i) from scale adjustment
		cosThetaIp = abs(cosThetaIp);
		fsum += Mp(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, bsdf.v[p]) * ap[p] * Np(phi, p, bsdf.s, bsdf.gammaO, gammaT);
	}

	// Compute contribution of remaining terms after P_MAX
	fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, bsdf.v[P_MAX]) * ap[P_MAX] * M_1_2PI;

	// wi is expressed the tangent plane of the hair
	// N == Z == vec3(0, 0, 1), cos(theta_i) == wi.z
	if (abs(wi.z) > 0) fsum /= abs(wi.z);

	return fsum;
}

vec3 BrdfHair::evalHair(const vec3& V,
	const vec3 L[3],
	int p,
	float& pdf) const
{
	vec4 baseColorBuffer(sigma_a.x, sigma_a.y, sigma_a.z, 1.0f);
	vec4 normalBuffer(1.0f, 0.0f, 0.0f, eta);
	vec4 materialBuffer(beta_m, beta_n, 0.0f, alpha);

	PrincipledHairBsdf bsdf;
	PrincipledHairBsdfDecodeGBuffer(h, baseColorBuffer, normalBuffer, materialBuffer, bsdf);

	Assert(bsdf.h == h);
	Assert(bsdf.eta == eta);
	Assert(bsdf.beta_m == beta_m);
	Assert(bsdf.beta_n == beta_n);
	Assert(bsdf.alpha == alpha);
	Assert(all(equal(bsdf.sigma_a, sigma_a)));

	// This is not correct, we should evaluate a specific part of hair bsdf based on p
	vec3 result = vec3(0.0f);
	pdf = 0.0f;

	return result;
}

void BrdfHair::sampleHair(const vec3& V,
	const float U1, const float U2,
	vec3 L[3]) const
{
	// TODO
}