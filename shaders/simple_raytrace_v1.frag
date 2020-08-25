#version 300 es

precision highp float;

uniform vec2 winsize;
uniform mat4 viewMatrix;
uniform mat4 projectionMatrix;
uniform vec3 cam_pos;

uniform float ratio;
uniform float scale;

#define EPSILON 1e-6
#define PI      3.14159265
#define PI_2    1.57079632

uniform int NB_SAMPLES;
uniform int NB_BOUNCES;
//uniform int NB_LIGHT_SAMPLES;

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

vec3 randHemiVector(vec3 normal)
{
    /*
    x = i / ( max_z * max_y )
    y = ( i / max_z ) % max_y
    z = i % max_z
*/
    int z = int( count / (uint(32) * uint(winsize.y))   );
    int y = int( (count / uint(32)) % uint(winsize.y)   );
    int x = int( count % uint(32)                       );
    //vec3 d = texelFetch(random_map, ivec3(x, y, z), 0).rgb;
    vec3 d = texelFetch(random_map, ivec3(gl_FragCoord.x, gl_FragCoord.y, count), 0).rgb;

    //count = (count+=uint(1)) %  uint(uint(winsize.x) * uint(winsize.y) * uint(32));
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

float angle(vec3 u, vec3 v)
{
    return dot(u,v)/(length(u)*length(v));
}

vec3 getPerpendicular(vec3 v)
{
    vec3 u = vec3(1);
    u.z = (v.x + v.y)/v.z;
    return u;
}

// ray origin, ray direction, sphere position, sphere radius
bool intersectSphere(vec3 Ro, vec3 Rd, vec3 Sp, float r, out Intersection inter)
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
                float t = t0;
                t0 = t1;
                t1 = t;
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
bool intersectPlane(vec3 Ro, vec3 Rd, vec3 p0, vec3 n, out Intersection inter)
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

// minimum and maximum coordinates of the cube
//                    o----X max
//                   /|   /|
//   y  z           o-o--o-o
//   | /            |/   |/
//   |/        min  X----o
//   o ---> x
//
// /!\ Normal at intersection is not computed yet !
//
bool intersectCube(vec3 Ro, vec3 Rd, vec3 min, vec3 max, out Intersection inter)
{
    float t0x = (min.x - Ro.x) / Rd.x; 
    float t1x = (max.x - Ro.x) / Rd.x; 
    float t0y = (min.y - Ro.y) / Rd.y; 
    float t1y = (max.y - Ro.y) / Rd.y; 
    float t0z = (min.z - Ro.z) / Rd.z; 
    float t1z = (max.z - Ro.z) / Rd.z;

    float tmin = t0x;
    float tmax = t1x;

    if(tmin > tmax) { swap(tmin, tmax); }
    if(t0y > t1y)   { swap(t0y, t1y); }

    if((tmin > t1y) || (t0y > tmax)) { return false; }


    if (t0y > tmin) 
        tmin = t0y; 
 
    if (t1y < tmax) 
        tmax = t1y; 
 
    if (t0z > t1z) { swap(t0z,t1z); } 
 
    if ((tmin > t1z) || (t0z > tmax)) 
        return false; 
    
    if(tmax < 0.0) // if intersection is behind
    {return false;}
    
    if (t0z > tmin) 
        tmin = t0z; 
 
    if (t1z < tmax) 
        tmax = t1z; 
    
    if(tmin > tmax)
        swap(tmax,tmin);

    inter.pos = Ro + Rd * tmin;
    inter.normal = vec3(0); ///////////////////////////!!!
    inter.t = tmin;
    return true; 
}

//
// three points describing the vector
//        
//          o v2
//         /|
//        / |
//       /  | 
//  v0  o---o v1
//
// use of inside/outside test for polygons
//
bool intersectTriangle(vec3 Ro, vec3 Rd, vec3 v0, vec3 v1, vec3 v2, out Intersection inter)
{
    // triangle's normal
    vec3 n  = cross(v1-v0,v2-v0);

    // d as in [ ax + by + cz + d = 0 ] 
    float d = dot(n,v0);

    // ray parallel to triangle ?
    float nRd = dot(n,Rd);
    if(nRd < EPSILON)
    {
        return false; 
    }

    // "distance" to the intersection from the ray's origin
    float t = - (dot(n, Ro) + d) / dot(n, Rd); 
    
    // triangle is behind the ray ?
    if(t < 0.0)
    {
        return false;
    }
    
    vec3 P = Ro + Rd * t;

    vec3 C;

    vec3 edge0 = v1 - v0;
    vec3 v0P = P - v0;
    C = cross(edge0,v0P);
    if(dot(n,C) < 0.0) { return false; }

    vec3 edge1 = v2 - v1;
    vec3 v1P = P - v1;
    C = cross(edge1,v1P);
    if(dot(n,C) < 0.0) { return false; }

    vec3 edge2 = v0 - v2;
    vec3 v2P = P - v2;
    C = cross(edge2,v2P);
    if(dot(n,C) < 0.0) { return false; }

    inter.pos = P;
    inter.normal = n;
    inter.t = t;

    return true;
}

bool trace(vec3 Ro, vec3 Rd, out Intersection closest)
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
    if(intersectPlane(Ro, Rd, vec3(0,0,0), vec3(0,1,0), inter))
    {
        if(inter.t < closest.t)
        {
            closest = inter;
            closest.matid = 2;
        }
        hit_smth = true;
    }
    return hit_smth;
}

float NormalDistribution(vec3 N, vec3 H, float alpha)
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

float GeometrySchlickGGX(float NdotV, float k)
{
    return NdotV / (NdotV * (1.0 - k) + k);
}

float Geometry(vec3 N, vec3 V, vec3 L, float alpha)
{
    float k = sqrt((2.0*alpha*alpha)/PI);
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx1 = GeometrySchlickGGX(NdotV, k);
    float ggx2 = GeometrySchlickGGX(NdotL, k);
	
    return ggx1 * ggx2;
}

// Using Fresnel-Schlick approximation
vec3 Fresnel(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

vec3 cookTorance(int lightIndex, vec3 Rd, Intersection inter)
{
    vec3 N = normalize(inter.normal);
    vec3 L = normalize(lights_pos[lightIndex] - inter.pos);
    vec3 V = normalize(camPos - inter.pos);
    vec3 H = normalize(V + L);

    vec3 albedo = mat_albedo[inter.matid];
    float alpha = mat_roughness[inter.matid];// * mat_roughness[inter.matid];
    float metallic = mat_metallic[inter.matid];

    float ND = NormalDistribution(N, H, alpha);
    float G = Geometry(N, V, L, alpha);
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
    //return PI_2 * specular;
}

// compute a rebound ray according to the view direction and the material[inter.matid]'s settings 
// S: random sample vector
vec3 sampleCookTorance(vec3 Sp, Intersection inter)
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
}

float orenNayar(int lightIndex, vec3 Rd, Intersection inter)
{
    float sigma = 2.0; // in radians
    vec3 N = normalize(inter.normal);
    vec3 L = normalize(lights_pos[lightIndex] - inter.pos);
    vec3 V = normalize(Rd);

    float NdotL = dot(L, N);

    if(NdotL < EPSILON)
    {
        return 0.0;
    }

    float LdotV = dot(L, V);
    float NdotV = dot(N, V);
    float albedo = 1.0;

    float s = LdotV - NdotL * NdotV;
    float t = mix(1.0, max(NdotL, NdotV), step(0.0, s));

    float sigma2 = sigma * sigma;
    float A = 1.0 + sigma2 * (albedo / (sigma2 + 0.13) + 0.5 / (sigma2 + 0.33));
    float B = 0.45 * sigma2 / (sigma2 + 0.09);

    return albedo * NdotL * (A + B * s / t) / PI;
}

float lambert(int lightIndex, vec3 Rd, Intersection inter)
{
    vec3 L = normalize(lights_pos[lightIndex] - inter.pos);
    vec3 N = normalize(inter.normal);
    
    float dist = length(L);
    float fade = 1.0/(2.0*dist*dist);
    fade = 1.0;
    return max( dot(N,L), 0.0)*fade;
}

vec3 samplePointLight(int lightIndex, vec3 Rd, in Intersection inter)
{
    vec3 L = lights_pos[lightIndex] - inter.pos; // to light center
    vec3 col;
    // if light has a size
    if(lights_size[lightIndex] > 0.0)
    {   
        vec3 sp = lights_pos[lightIndex] + (normalize(randHemiVector(-L)) * lights_size[lightIndex]);
        vec3 Lsp = sp - inter.pos;
        float radius = lights_size[lightIndex];

        Intersection itt;
        if(!trace(inter.pos, Lsp, itt))
        {
            float distance = length(Lsp);
            float attenuation = 4.0 * PI * distance * distance;
            vec3 radiance = lights_col[lightIndex] * (lights_intensity[lightIndex] / attenuation);
            col += radiance;
        }

        return col;
    }

    else
    {
        Intersection itt;
        if(!trace(inter.pos, L, itt))
        {
            float distance = length(L);
            float attenuation = 1.0 / (4.0 * PI * distance * distance);
            vec3 radiance = lights_col[lightIndex] * attenuation * lights_intensity[lightIndex];

            return radiance;;
        }
    }

    return vec3(0);
}

void main()
{
    //count = uint(gl_FragCoord.x) + uint(gl_FragCoord.y) * uint(winsize.x);
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
                col += cookTorance(i, Rd, interPrimary) * samplePointLight(i, Rd, interPrimary);
            }

            // SECONDARY RAYS

            // the intersection from which we are looking
            Intersection interPrev;
            // the intersection we are looking at with Rds
            Intersection interNext;

            interPrev = interPrimary;
            
            vec3 acc;// = vec3(mat_albedo[interPrev.matid]);
            
            for(int j = 0; j < NB_BOUNCES ; j++)
            {
                // compute secondary ray
                // random ray from interPrev.position
                vec3 Rds = normalize(randHemiVector(interPrev.normal));
                float cosTheta = dot(interPrev.normal, Rds);
                vec3 acc_current = vec3(0);

                // intersection test
                // if we hit an object
                if(trace(interPrev.pos + normalize(interPrev.normal)*EPSILON, Rds, interNext))
                {
                    // compute direct illumination
                    for(int i = 0; i < NB_LIGHTS ; i++)
                    {
                        acc_current += cookTorance(i, Rd, interNext) * samplePointLight(i, Rds, interNext);
                    }
                }
                else // nothing is hit
                {
                    if(ambiant_contribute) 
                    {
                        acc_current += ambiant_light * ambiant_intensity;
                    }
                    // stopping the loop
                    j = NB_BOUNCES;
                }

                // accumulate indirect lighting
                acc += acc_current * cosTheta;

                //update interSecondary
                interPrev = interNext;

            } // end for bounces
            
            col += acc;
            col *= mat_albedo[interPrimary.matid] / PI;

        } // end if(hit primary)

        // if primary intersection fails
        else{
            if(ambiant_contribute) col += ambiant_light * ambiant_intensity;
        }
        
        // accumulate light from paths
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