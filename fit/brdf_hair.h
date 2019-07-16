#ifndef _BRDF_HAIR_
#define _BRDF_HAIR_

#include "brdf.h"

class BrdfHair : public Brdf
{
public:

	virtual float eval(const vec3& V, const vec3& L, const float alpha, float& pdf) const
	{
		// Original Brdf::eval
		return 0.0f;
	}

    virtual vec3 sample(const vec3& V, const float alpha, const float U1, const float U2) const
    {
		// Original Brdf::sample
		return vec3(0.0f);
    }

	// Evaluate hair bsdf to get an idea of direction distribution
	// See https://www.pbrt.org/hair.pdf.
	// Hair bsdf is factored into 3 terms: Mp, Ap and Np. Don't know if I could remove some,
	virtual vec3 evalHair(const vec3& V,
		const vec3& L,
		const vec3& Tangent, // this vector probably can be removed as we can assume we will be in BRDF space and X could be our default orientation
		float h,
		const vec3& sigma_a,
		float eta,
		float beta_m,
		float beta_n,
		float alpha,
		int p, // 0 == R, 1 == TT, 2 == TRT, 3 == whatever is left
		float& pdf) const;

	// Get importance sampled lobes to fit
	virtual void sampleHair(const vec3& V,
		const float U1, const float U2,
		const vec3& Tangent, // this vector probably can be removed as we can assume we will be in BRDF space and X could be our default orientation
		float h, // hair parameters, don't know if we could remove some
		const vec3& sigma_a,
		float eta,
		float beta_m,
		float beta_n,
		float alpha,
		vec3 L[3]) const; // Get R, TT and TRT vectors, each of them will get fitted individually
};

#endif
