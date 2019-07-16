#ifndef _BRDF_
#define _BRDF_

#include <glm/glm.hpp>
using namespace glm;

class Brdf
{
public:
    // evaluation of the cosine-weighted BRDF
    // pdf is set to the PDF of sampling L
    virtual float eval(const vec3& V, const vec3& L, const float alpha, float& pdf) const = 0;

    // sampling
    virtual vec3 sample(const vec3& V, const float alpha, const float U1, const float U2) const = 0;


	// EXPERIMENTATION: SEE HAIRBSDF.H
	virtual vec3 evalHair(const vec3& V,
		const vec3& L,
		const vec3& Tangent,
		float h,
		const vec3& sigma_a,
		float eta,
		float beta_m,
		float beta_n,
		float alpha,
		int p,
		float& pdf) const { return vec3(0.0f); }

	virtual void sampleHair(const vec3& V,
		const float U1, const float U2,
		const vec3& Tangent,
		float h,
		const vec3& sigma_a,
		float eta,
		float beta_m,
		float beta_n,
		float alpha,
		vec3 L[3]) const {}
};

#endif
