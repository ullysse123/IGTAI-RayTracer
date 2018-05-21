
#include "raytracer.h"
#include "scene_types.h"
#include "ray.h"
#include "image.h"
#include "kdtree.h"
#include <stdio.h>

#define PI 3.141592653f
#define MAX_DEPTH 1

/// acne_eps is a small constant used to prevent acne when computing intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;

bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj) {
  //!\todo
  bool intersect;
  float t;
  vec3 dir = ray->dir;
  vec3 n = obj->geom.plane.normal;
  point3 O = ray->orig;
  
  if ( dot(dir,n) == 0 ){
      intersect = false;
  }else{
      intersect = true;
      t = -( dot(O,n) + obj->geom.plane.dist)/(dot(dir,n));
      if ( !(ray->tmin < t && ray->tmax > t) ){
          intersect = false;
      }
  }
  
  if (intersect){
      intersection->mat = &obj->mat;
      intersection->position = rayAt(*ray,t);
      intersection->normal = n;
      ray->tmax = t;
  }
  
  return intersect;

}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj) {
  //!\todo
  //t² + 2t(vec(d).(O-C))+((O-C).(O-C)-R²)
  bool intersect;
  
  point3 O = ray->orig;
  point3 C = obj->geom.sphere.center;
  float radius = obj->geom.sphere.radius;
  vec3 dir = ray->dir;
  float t;
  
  float a = 1.0f;
  float b = 2.0f*(dot(dir,(O-C)));
  float c = dot((O-C),(O-C))-radius*radius;
  
  float delta = b*b - 4.0f * a * c;
  
  if (delta > 0){
    float tUn = (-b+sqrt(delta))/(2.0f*a);
    float tDeux = (-b-sqrt(delta))/(2.0f*a);
    bool un = (tUn >= ray->tmin && tUn <= ray->tmax);
    bool deux = (tDeux >= ray->tmin && tDeux <= ray->tmax);
    if(un && deux){
      if(tUn>tDeux){
        t = tDeux;
        intersect = true;
      }else{
        t = tUn;
        intersect = true;
      }
    }else{
      if(un && !deux){
        t = tUn;
        intersect = true;
      }else{
        if(!un && deux){
          t = tDeux;
          intersect = true;
        }else{
          intersect = false;
        }
      }
    }
  }else{
    if(delta == 0){
      //Une solution
      t = -b/(2.0f*a);
      intersect = (t >= ray->tmin && t <= ray->tmax);
      
    }else{
      //Aucune solution
      intersect = false;
    }
  }
  
  if(intersect){
    intersection->mat = &obj->mat;
    intersection->position = rayAt(*ray,t);
    intersection->normal = normalize(intersection->position - C);
    ray->tmax = t;
  }
  
  return intersect;
}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
  bool hasIntersection = false;
  size_t objectCount = scene->objects.size();
  
  //!\todo loop on each object of the scene to compute intersection
  for(size_t i=0; i<objectCount; i++){
      switch(scene->objects[i]->geom.type){
          case SPHERE : 
            hasIntersection = (intersectSphere(ray,intersection,scene->objects[i])) || hasIntersection;
              break;
          case PLANE : 
            hasIntersection = (intersectPlane(ray,intersection,scene->objects[i])) || hasIntersection;
              break;
      }
  }
  return hasIntersection;
}

/* --------------------------------------------------------------------------- */
/*
 *	The following functions are coded from Cook-Torrance bsdf model description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba renderer)
 */

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */
float RDM_Beckmann(float NdotH, float alpha) {

  //! \todo compute Beckmann normal distribution
  float cos2 = NdotH * NdotH;
  float numerateur = exp(-((1-cos2)/cos2)/(alpha*alpha));
  float denumerateur = PI*alpha*alpha*cos2*cos2;
  return numerateur/denumerateur;

}

// Fresnel term computation. Implantation of the exact computation. we can use the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {

  //! \todo compute Fresnel term
  float cos2 = LdotH * LdotH;
  float sin2 = ((extIOR/intIOR)*(extIOR/intIOR))*(1-cos2);
  if (sin2 > 1) return 1;
  
  float rs = ((extIOR*LdotH - intIOR*sqrt(1 - sin2))*(extIOR*LdotH - intIOR*sqrt(1 - sin2)))/((extIOR*LdotH + intIOR*sqrt(1-sin2))*(extIOR*LdotH + intIOR*sqrt(1-sin2)));
  float rp = ((extIOR*sqrt(1 - sin2)-intIOR*LdotH)*(extIOR*sqrt(1 - sin2)-intIOR*LdotH))/((extIOR*sqrt(1 - sin2)+intIOR*LdotH)*(extIOR*sqrt(1 - sin2)+intIOR*LdotH));
  
  return (rs + rp)/2;

}


// Shadowing and masking function. Linked with the NDF. Here, Smith function, suitable for Beckmann NDF
float RDM_chiplus(float c) {
  return (c > 0.f) ? 1.f : 0.f;
}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha) {

  //!\todo compute G1 term of the Smith fonction
  float ret = 1.f;
  float k = DdotH / DdotN;
  float cos2 = DdotN*DdotN;
  float b = 1 /(alpha*(sqrt(1-cos2)/DdotN));
  
  if(k>0.f && b<1.6f){
    ret = (3.535f*b + 2.181f*b*b)/(1.0f + 2.276f*b + 2.577f*b*b);
  }
  
  ret = ret * RDM_chiplus(k);
  
  return ret;

}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN, float alpha) {

  //!\todo the Smith fonction
  return RDM_G1(LdotH,LdotN,alpha) * RDM_G1(VdotH,VdotN,alpha);


}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m) {

  //!\todo specular term of the bsdf, using D = RDB_Beckmann, F = RDM_Fresnel, G = RDM_Smith
  float D = RDM_Beckmann(NdotH,m->roughness);
  float F = RDM_Fresnel(LdotH,1,m->IOR);
  float G = RDM_Smith(LdotH,LdotN,VdotH,VdotN,m->roughness);
  return (m->specularColor*D*F*G)/(4*LdotN*VdotH);

  
}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {

  //!\todo compute diffuse component of the bsdf
  return m->diffuseColor/PI;

}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m) {

  //! \todo compute bsdf diffuse and specular term
  return RDM_bsdf_d(m)+RDM_bsdf_s(LdotH,NdotH,VdotH,LdotN,VdotN,m);

}




/* --------------------------------------------------------------------------- */

color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat ){
  color3 ret = color3(0.f);

  //! \todo compute bsdf, return the shaded color taking into account the
  //! lightcolor
  if(dot(l,n)>=0.f){
    //ret = (mat->diffuseColor/PI)*dot(l,n)*lc;
    vec3 h = normalize((v+l)/length(v+l));
    ret = RDM_bsdf(dot(l,h),dot(n,h),dot(v,h),dot(l,n),dot(v,n),mat)*lc*dot(l,n);
  }

  return ret;
	    
}

//! if tree is not null, use intersectKdTree to compute the intersection instead of intersect scene
color3 trace_ray(Scene * scene, Ray *ray, KdTree *tree) {
  color3 ret = color3(0.f,0.f,0.f);
  
  if(ray->depth > MAX_DEPTH) return ret;
  
  Intersection intersection;
  Ray r;
  if(!intersectScene(scene,ray,&intersection)){
    
    ret = scene->skyColor;
    
  }else{
    
    //ret = intersection.normal*0.5f+0.5f;
    for(Light *light : scene->lights){
        vec3 direction_lum = light->position - intersection.position;
        vec3 l = normalize(direction_lum);
        rayInit(&r,intersection.position + acne_eps*l,l, acne_eps, length(direction_lum));
        Intersection ombre;
        if(!intersectScene(scene,&r,&ombre)){
            ret = ret + shade(intersection.normal, -(ray->dir), l, light->color, intersection.mat);
        }
    }
    vec3 refDir = normalize(reflect(ray->dir,intersection.normal));
    float LdotH = dot(refDir,normalize(-ray->dir + refDir));
    rayInit(ray,intersection.position,refDir,acne_eps,100000,ray->depth+1);
    ret = ret + RDM_Fresnel(LdotH,1,intersection.mat->IOR) * trace_ray(scene,ray,tree);
    
  }
  
  return ret;
}

void renderImage(Image *img, Scene *scene) {

  //! This function is already operational, you might modify it for antialiasing and kdtree initializaion
  float aspect = 1.f/scene->cam.aspect;
    
  KdTree *tree =  NULL;


  //! \todo initialize KdTree

  float delta_y = 1.f / (img->height * 0.5f); //! one pixel size
  vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step 
  vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) * aspect * scene->cam.ydir;

  float delta_x = 1.f / (img->width * 0.5f);
  vec3 dx = delta_x * scene->cam.xdir;
  vec3 ray_delta_x = (0.5f - img->width * 0.5f) / (img->width * 0.5f) *scene->cam.xdir;
  
    
  for(size_t j=0; j<img->height; j++) {
    if(j!=0) printf("\033[A\r");
    float progress = (float)j/img->height*100.f;
    printf("progress\t[");
    int cpt = 0;
    for(cpt = 0; cpt<progress; cpt+=5) printf(".");
    for(       ; cpt<100; cpt+=5) printf(" ");
    printf("]\n");
#pragma omp parallel for
    for(size_t i=0; i<img->width; i++) {
      color3 *ptr = getPixelPtr(img, i,j);
      vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y + float(i)*dx + float(j)*dy;

      Ray rx;
      rayInit(&rx, scene->cam.position, normalize(ray_dir));
      *ptr = trace_ray(scene, &rx, tree);

    }
  }
}
