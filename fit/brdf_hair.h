#ifndef _BRDF_HAIR_
#define _BRDF_HAIR_

#include "brdf.h"

class BrdfHair : public Brdf
{
public:

    virtual float eval(const vec3& V, const vec3& L, const float alpha, float& pdf) const
    {
		return 0.0f;
    }

    virtual vec3 sample(const vec3& V, const float alpha, const float U1, const float U2) const
    {
		return vec3(0.0f);
    }
};

#endif
