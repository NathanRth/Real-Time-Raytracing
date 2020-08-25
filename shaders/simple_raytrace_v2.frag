#version 300 es

precision highp float;

uniform vec2 winsize;
uniform mat4 viewMatrix;
uniform mat4 projectionMatrix;
uniform vec3 cam_pos;

uniform float ratio;
uniform float scale;

#define EPSILON 1e-6
#define PI 3.14159265

uniform int NB_SAMPLES;
uniform int NB_BOUNCES;
uniform int NB_LIGHT_SAMPLES;

#define RAND_MAP_SIZE   2000
uniform highp sampler3D random_map;
uniform bool show_random;

uniform vec3 ambiant_light;
uniform bool ambiant_contribute;
uniform float ambiant_intensity;

#define MAX_NB_LIGHTS 2
uniform int NB_LIGHTS;
uniform vec3 lights_pos[MAX_NB_LIGHTS];
uniform vec3 lights_col[MAX_NB_LIGHTS];
uniform float lights_intensity[MAX_NB_LIGHTS];
uniform float lights_size[MAX_NB_LIGHTS];

#define MAX_NB_SPHERES 3
uniform int NB_SPHERES;
uniform vec3 spheres_pos[MAX_NB_SPHERES];
uniform float spheres_radius[MAX_NB_SPHERES];
uniform int spheres_matid[MAX_NB_SPHERES];

#define MAX_NB_PLANES 3
uniform int NB_PLANES;
uniform vec3 planes_pos[MAX_NB_PLANES];
uniform vec3 planes_normal[MAX_NB_PLANES];
uniform int planes_matid[MAX_NB_PLANES];

#define MAX_NB_MAT   4
uniform int NB_MATERIALS;
uniform vec3 mat_albedo[MAX_NB_MAT];
uniform float mat_metallic[MAX_NB_MAT];
uniform float mat_roughness[MAX_NB_MAT];

vec3 camPos;

struct Intersection
{
    vec3 pos;       // 3D position of intersection
    vec3 normal;    // normal of the surface at intersection
    float t;        // "distance" of the intersection from Ro
    int matid;      // material id of intersected surface
};

out vec4 frag_out;
uint count;

vec3 randHemiVector(in vec3 normal)
{
    vec3 d = texelFetch(random_map, ivec3(gl_FragCoord.x, gl_FragCoord.y, count), 0).rgb;
    count = (count+=uint(1)) %  uint(32);
    
    // if d is not on the same "side" (same hemisphere oriented along <normal>), flip d
    if(dot(normal,d)>EPSILON)
    {
        return d;
    }
    else
    {
        return -d;
    }
}

void swap(inout float a, inout float b)
{
    float t = a;
    a = b;
    b = t;
}

// ray origin, ray direction, sphere position, sphere radius
bool intersectSphere(in vec3 Ro, in vec3 Rd, in vec3 Sp, in float r, out Intersection inter)
{
    vec3 L = Ro - Sp;
    
    float a = dot(Rd,Rd);
    float b = 2.0 * dot(Rd,L);
    float c = dot(L,L)-r*r;

    float delta = b*b - 4.0*a*c;

    if (delta < 0.0)
    {
        return false;
    }

    if(delta < 1e-6)
    {
        float t = - b / (2.0 * a);
        if(t > EPSILON) // if intersection is behind
        {
            inter.pos = Ro + Rd * t;
            inter.normal = normalize(inter.pos - Sp);
            inter.t = t;
            return true;
        }
    }
    else
    {
        float t0 = - (b + sqrt(delta)) / (2.0 * a);
        float t1 = - (b - sqrt(delta)) / (2.0 * a);
        if(t0 > t1)
        {
            swap(t0,t1);
        }
        if(t0 > EPSILON) // if intersection is behind
        {
            if(t0 > t1)
            {
                swap(t0,t1);
            }
            inter.pos = Ro + Rd * t0;
            inter.normal = normalize(inter.pos - Sp);
            inter.t = t0;
            return true;
        }
    }
    return false;
}

// ray origin, ray direction, point of plane, plane normal
bool intersectPlane(in vec3 Ro, in vec3 Rd, in vec3 p0, in vec3 n, out Intersection inter)
{
    float nRd = dot(n,Rd);
    if(abs(nRd) > EPSILON )
    {
        vec3 p0Ro = p0 - Ro;
        float t = dot(p0Ro, n) / nRd;
        if( t > EPSILON) // if intersection is behind
        {
            inter.normal = n;
            inter.pos = Ro + Rd * t;
            inter.t = t;
            return true;
        }
    }
    return false;
}

bool trace(in vec3 Ro, in vec3 Rd, out Intersection closest)
{
    Intersection inter;
    closest.t = 10000.0;
    bool hit_smth = false;

    for(int i = 0; i < NB_SPHERES ; i++)
    {
        if(intersectSphere(Ro, Rd,spheres_pos[i], spheres_radius[i], inter))
        {
            if(inter.t < closest.t)
            {
                closest = inter;
                closest.matid = spheres_matid[i];
            }
            hit_smth = true;
        }
    }
    for(int i = 0; i < NB_PLANES ; i++)
    {
        if(intersectPlane(Ro, Rd, planes_pos[i], planes_normal[i], inter))
        {
            if(inter.t < closest.t)
            {
                closest = inter;
                closest.matid = planes_matid[i];
            }
            hit_smth = true;
        }
    }
    return hit_smth;
}

float NormalDistribution(in vec3 N, in vec3 H, in float alpha)
{
    float alpha2 = alpha * alpha;
    float NdotH = max(dot(N,H), 0.0);
    float NdotH2 = NdotH * NdotH;
    float denom = NdotH2 * (alpha2 - 1.0)+1.0;
    return alpha2 / (PI*denom*denom);
/*
    float NdotH = dot(N,H);
    float NdotH2 = NdotH * NdotH;
    float NdotH4 = NdotH2 * NdotH2;
    return (1.0 / PI * alpha2 * NdotH4) * exp((NdotH2-1.0)/alpha2*NdotH2);*/
}

float GeometrySchlickGGX(in float NdotV, in float k)
{
    return NdotV / (NdotV * (1.0 - k) + k);
}

float Geometry(in vec3 N, in vec3 V, in vec3 L, in float k)
{
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx1 = GeometrySchlickGGX(NdotV, k);
    float ggx2 = GeometrySchlickGGX(NdotL, k);
	
    return ggx1 * ggx2;
}

// Using Fresnel-Schlick approximation
vec3 Fresnel(in float cosTheta, in vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}


// Evaluation of the CookTorance BRDF at inter.pos according to lightPos
vec3 cookTorance(in vec3 lightPos, in Intersection inter)
{
    // normal at inter.pos
    vec3 N = normalize(inter.normal);
    // lighting vector (toward the light)
    vec3 L = normalize(lightPos - inter.pos);
    // viewing vector (towards the viewer)
    vec3 V = normalize(camPos - inter.pos);
    // halway vector
    vec3 H = normalize(V + L);

    vec3 albedo = mat_albedo[inter.matid];
    float roughness = mat_roughness[inter.matid];
    float metallic = mat_metallic[inter.matid];

    float ND = NormalDistribution(N, H, roughness);
    float G = Geometry(N, V, L, roughness);
    vec3 F0 = vec3(0.04);
    F0 = mix(F0, albedo, metallic);
    vec3 F = Fresnel(max(dot(H,V),0.0),F0);

    vec3 num = ND * F * G;
    float denom = 4.0 * max(dot(N,V), 0.0) * max(dot(N, L), 0.0);
    vec3 specular = num / max(denom,EPSILON);

    vec3 ks = F;
    vec3 kd = vec3(1.0)-ks;
    kd *= 1.0 - metallic;

    float NdotL = max(dot(N,L), 0.0);
    return (kd * albedo / PI + specular) * NdotL;
}

/*vec3 pdfCookTorance(in vec3 Sp, in Intersection inter)
{
    vec3 N = normalize(inter.normal);
    vec3 S = normalize(Sp);
    vec3 V = normalize(camPos - inter.pos);
    vec3 H = normalize(V + S);

    float roughness = mat_roughness[inter.matid];
    float metallic = mat_metallic[inter.matid];
    vec3 albedo = mat_albedo[inter.matid];

    float ND = NormalDistribution(N, H, roughness);
    float G = Geometry(N, V, S, roughness);
    vec3 F0 = mix(vec3(0.04),albedo,metallic);
    vec3 F = Fresnel(max(dot(H,V),0.0),F0);

    return (PI / 2.0)*(ND*F*G/dot(N,V));
}*/

vec3 samplePointLight(in int lightIndex, vec3 Rd, in Intersection inter)
{
    vec3 N = inter.normal;
    vec3 L = lights_pos[lightIndex] - inter.pos; // to light center
    vec3 col;
    // if light has a size
    if(lights_size[lightIndex] > 0.0)
    {
        
        for(int i = 0 ; i < NB_LIGHT_SAMPLES ; i++)
        {    
            vec3 sp = lights_pos[lightIndex] + (normalize(randHemiVector(-L)) * lights_size[lightIndex]);
            vec3 Lsp = sp - inter.pos;
            float NdotLsp = dot(normalize(N),normalize(Lsp));

            Intersection itt;
            if(!trace(inter.pos, Lsp, itt))
            {
                float distance = length(Lsp);
                float attenuation = 1.0 / (distance * distance);
                vec3 radiance = lights_col[lightIndex] * attenuation * lights_intensity[lightIndex];// * NdotLsp;
                col += radiance;
            }
        }
        return col/float(NB_LIGHT_SAMPLES);
    }

    else
    {
        L = lights_pos[lightIndex] - inter.pos;

        Intersection itt;
        if(!trace(inter.pos, L, itt))
        {
            float distance = length(L);
            float attenuation = 1.0 / (distance * distance);
            vec3 radiance = lights_col[lightIndex] * attenuation * lights_intensity[lightIndex];

            return radiance;
        }
    }

    return vec3(0);
}

void main()
{
    //count = uint(gl_FragCoord.x) + uint(gl_FragCoord.y) * uint(winsize.x) * uint(128);
    count = uint(0);
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;
    
    float Px = (2.0 * x) / winsize.x -1.0;
    float Py = -(1.0 - (2.0 * y) / winsize.y);
    
    vec4 ray_clip = vec4(Px, Py, -1.0, 1.0);
    vec4 ray_eye = inverse(projectionMatrix) * ray_clip;
    ray_eye = vec4(ray_eye.xy, -1.0, 0.0);
    vec3 ray_world = (inverse(viewMatrix) * ray_eye).xyz;
    
    vec4 rayOrigin = inverse(viewMatrix)*vec4(0,0,0,1.0);
    camPos = rayOrigin.xyz;

    vec3 Rd = normalize(ray_world);
    vec3 Ro = camPos;

    for(int s = 0 ; s < NB_SAMPLES ; s++)
    {
        vec3 col = vec3(0);

        Intersection interPrimary;
        
        if(trace(Ro, Rd, interPrimary))
        {
            for(int i = 0; i < NB_LIGHTS ; i++)
            {
                vec3 radiance = samplePointLight(i, Rd, interPrimary);
                col += cookTorance(lights_pos[i], interPrimary) * radiance;
                //col += radiance;
            }

            // SECONDARY RAYS

            Intersection interPrev;
            Intersection interNext;

            interPrev = interPrimary;
            
            //vec3 acc = vec3(mat_albedo[interPrev.matid]);
            
            for(int j = 0; j < NB_BOUNCES ; j++)
            {
                // compute secondary ray
                vec3 Rds = randHemiVector(interPrev.normal);
                
                //vec3 acc_current = vec3(0);

                //send secondary ray
                if(trace(interPrev.pos + normalize(interPrev.normal)*EPSILON, Rds, interNext))
                {
                    for(int i = 0; i < NB_LIGHTS ; i++)
                    {
                        // for each light at secondary intersection
                    }
                }
                else // nothing is hit
                {
                    if(ambiant_contribute) 
                    {
                        //acc_current += ambiant_light * ambiant_intensity;
                    }
                    //acc *= acc_current;
                    break;
                }

                //acc *= acc_current;

                //update interSecondary
                interPrev = interNext;

            } // end for bounces
            
            //col += acc;
            //col *= mat_albedo[interPrimary.matid];
            
        } // end if(hit primary)
        else{
            if(ambiant_contribute) col += ambiant_light * ambiant_intensity;
        }
        
        frag_out += vec4(col,1.0);
    
    } // end for(samples)
    
    if(show_random)
    {
        frag_out = vec4(randHemiVector(vec3(1,0,0)),1.0);
    }
    else
    {
        frag_out /= float(NB_SAMPLES);
    }
}