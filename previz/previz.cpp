//////////////////////////////////////////////////////////////////////////////////
// This is a front end for a set of viewer clases for the Carnegie Mellon
// Motion Capture Database: 
//    
//    http://mocap.cs.cmu.edu/
//
// The original viewer code was downloaded from:
//
//   http://graphics.cs.cmu.edu/software/mocapPlayer.zip
//
// where it is credited to James McCann (Adobe), Jernej Barbic (USC),
// and Yili Zhao (USC). There are also comments in it that suggest
// and Alla Safonova (UPenn) and Kiran Bhat (ILM) also had a hand in writing it.
//
//////////////////////////////////////////////////////////////////////////////////
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <float.h>
#include "SETTINGS.h"
#include "skeleton.h"
#include "displaySkeleton.h"
#include "motion.h"

using namespace std;

// Stick-man classes
DisplaySkeleton displayer;    
Skeleton* skeleton;
Motion* motion;

// cylinder class

struct Cylinder {
  VEC3 center;
  VEC3 axis;

  VEC3 left;
  VEC3 right;

  const float height;
  const float radius;
};

struct Light {
  VEC3 position;
  VEC3 color;
};

struct Triangle {
  VEC3 a;
  VEC3 b;
  VEC3 c;

  bool shiny;
};

int windowWidth = 640;
int windowHeight = 480;

VEC3 eye(-5, 1, 5);

// VEC3 eye(-5, 0.5, -5);

VEC3 lookingAt(3, 0.5, 1);
VEC3 up(0,1,0);

// scene geometry
// spheres
vector<VEC3> sphereCenters;
vector<float> sphereRadii;
vector<VEC3> sphereColors;

// cylinders
vector<Cylinder> cylinders;

// lights
vector<Light> lights;

//triangles
vector<Triangle> triangles;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void writePPM(const string& filename, int& xRes, int& yRes, const float* values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
bool raySphereIntersect(const VEC3& center, 
                        const float radius, 
                        const VEC3& rayPos, 
                        const VEC3& rayDir,
                        float& t)
{
  const VEC3 op = center - rayPos;
  const float eps = 1e-8;
  const float b = op.dot(rayDir);
  float det = b * b - op.dot(op) + radius * radius;

  // determinant check
  if (det < 0)
    return false; 
  
  det = sqrt(det);
  t = b - det;
  if (t <= eps)
  {
    t = b + det;
    if (t <= eps)
      t = -1;
  }

  if (t < 0) return false;
  return true;
}

bool rayCylinderIntersect(const float radius, const VEC3& axis, const VEC3& center,
                        const VEC3& rayPos, const VEC3& left, const VEC3& right,
                        const VEC3& rayDir, float& t)
{
  const float a = rayDir.dot(rayDir) - pow(rayDir.dot(axis), 2);
  const float b = 2 * (rayDir.dot(rayPos - center) - (rayDir.dot(axis) * (rayPos - center).dot(axis)));
  const float c = (rayPos - center).dot(rayPos - center) - pow((rayPos-center).dot(axis), 2) - (radius * radius);

  float det = b * b - (4 * a * c);
  float yCheck;
  
  if (det < 0) return false;

  det = sqrt(det);

  t = (-b - det) / (2 * a);

  if (t < 0) return false;

  float max = sqrt(pow((left-right).norm() / 2, 2) + pow(radius, 2)); //pythagoras theorem
  VEC3 point = rayPos + rayDir * t;
  VEC3 len = point - center;

  if (len.norm() > max) // if t1 is too high we try with t2
    t = (-b + det) / (2 * a);

  point = rayPos + rayDir * t;
  len = point - center;

  if (len.norm() > max)
    return false;

  return true;
}

bool rayTriangleIntersect(Triangle tri, const VEC3& rayPos, const VEC3& rayDir, float& t) {
  VEC3 a = tri.b - tri.a; // edge 0
  VEC3 b = tri.c - tri.a; // edge 1
  VEC3 c = a.cross(b).normalized(); // this is the triangle's normal
  float NdotRayDirection = c.dot(rayDir);
  if (fabs(NdotRayDirection) < 0.0008)
    return false; // almost 0, parallel don't intersect

  float d = -c.dot(tri.a);
  t = -(c.dot(rayPos) + d) / c.dot(rayDir);

  if (t < 0) return false; // the triangle is behind

  VEC3 pointHit = rayPos + t * rayDir;
  VEC3 C;


  // edge 0
  VEC3 edge0 = tri.b - tri.a; 
  VEC3 vp0 = pointHit - tri.a;
  C = edge0.cross(vp0);
  if (c.dot(C) < 0) return false; // P is on the right side

  // edge 1
  VEC3 edge1 = tri.c - tri.b; 
  VEC3 vp1 = pointHit - tri.b;
  C = edge1.cross(vp1);
  if (c.dot(C) < 0)  return false; // P is on the right side

  // edge 2
  VEC3 edge2 = tri.a - tri.c; 
  VEC3 vp2 = pointHit - tri.c;
  C = edge2.cross(vp2);
  if (c.dot(C) < 0) return false; // P is on the right side;

  return true;
}
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

bool getInShadow(VEC3 point, VEC3 direction) {
  float tMin = FLT_MAX;
  for (int y = 0; y < sphereCenters.size(); y++)
  {
    if (raySphereIntersect(sphereCenters[y], sphereRadii[y], point, direction, tMin))
    {
      return true;
    }
  }
  for (int z = 0; z < cylinders.size(); z++)
  {
    Cylinder currentCylinder = cylinders[z];
    if (rayCylinderIntersect(currentCylinder.radius, currentCylinder.axis, currentCylinder.center, point, currentCylinder.left, currentCylinder.right, direction, tMin))
    {
      return true;
    }
  }
  for (int a = 0; a < triangles.size(); a++)
  {
    Triangle currentTriangle = triangles[a];
    if (rayTriangleIntersect(currentTriangle, point, direction, tMin))
    {
      return true;
    }
  }
  return false;
}

void rayColor(const VEC3& rayPos, const VEC3& rayDir, VEC3& pixelColor) 
{
  // pixelColor = VEC3(0.5294,0.8078,0.921); //sky blue
  pixelColor = VEC3(0,0,0);
  VEC3 c(0, 0, 0);
  VEC3 boneColor(1, 1, 1);

  // look for intersections
  bool sphereHit = false;
  bool cylinderHit = false;
  bool triangleHit = false;

  bool shadowFlag = false;

  bool shiny = false;

  int hitID = -1;
  float tMinFound = FLT_MAX;

  VEC3 point;
  VEC3 normal;

  for (int a = 0; a < triangles.size(); a++)
  {
    float tMin = FLT_MAX;
    if (rayTriangleIntersect(triangles[a], rayPos, rayDir, tMin))
    { 
      // is the closest so far?
      if (tMin < tMinFound)
      {
        triangleHit = true;
        sphereHit = false;
        cylinderHit = false;
        
        tMinFound = tMin; 
        point = (rayPos * .9995) + (rayDir * tMinFound);
        normal = (triangles[a].a).cross(triangles[a].b).normalized(); // this is the triangle's normal

        hitID = a;
      }
    }
  }

  for (int y = 0; y < sphereCenters.size(); y++)
  {
    float tMin = FLT_MAX;
    if (raySphereIntersect(sphereCenters[y], sphereRadii[y], rayPos, rayDir, tMin))
    { 
      // is the closest so far?
      if (tMin < tMinFound)
      {
        sphereHit = true;
        triangleHit = false;
        cylinderHit = false;

        tMinFound = tMin;

        point = (rayPos * 1.0008) + (rayDir * tMinFound);
        normal = (point - sphereCenters[y]) / sphereRadii[y];

        hitID = y;
      }
    }
  }

  for (int z = 0; z < cylinders.size(); z++)
  {
    float tMin = FLT_MAX;
    Cylinder currentCylinder = cylinders[z];
    if (rayCylinderIntersect(currentCylinder.radius, currentCylinder.axis, currentCylinder.center, rayPos, currentCylinder.left, currentCylinder.right, rayDir, tMin))
    { 
      // is the closest so far?
      if (tMin < tMinFound)
      {
        sphereHit = false;
        triangleHit = false;
        cylinderHit = true;
        tMinFound = tMin;

        point = (rayPos * 1.0008) + (rayDir * tMinFound);

        float m = rayDir.dot(currentCylinder.axis) * tMinFound + (rayPos - currentCylinder.center).dot(currentCylinder.axis);

        normal = (point - currentCylinder.left - currentCylinder.axis * m).normalized();

        hitID = z;
      }
    }
  }
  
  if (triangleHit) {
    for(int l = 0; l < lights.size(); l++) {
      VEC3 to_light = ((eye - point) + (lights[l].position - point)).normalized();

      float l_norm = normal.dot((lights[l].position - point).normalized());
      float h_norm = normal.dot(to_light);

      VEC3 highlight = ((-1 * (lights[l].position - point).normalized()) + 2 * l_norm * normal).normalized();
      VEC3 phong = lights[l].color * pow(max(float(0), float(highlight.dot(-rayDir))), 10);

      shadowFlag = !getInShadow(point, (lights[l].position - point).normalized());
      c += shadowFlag * (boneColor.cwiseProduct(max(float(0), l_norm) * lights[l].color + phong));    
    }
    pixelColor = c;
  }

  // No intersection, return white
  if (hitID == -1)
    return;

  if (sphereHit) {
    for(int l = 0; l < lights.size(); l++) {
      VEC3 to_light = ((eye - point) + (lights[l].position - point)).normalized();

      float l_norm = normal.dot((lights[l].position - point).normalized());
      float h_norm = normal.dot(to_light);

      VEC3 highlight = ((-1 * (lights[l].position - point).normalized()) + 2 * l_norm * normal).normalized();
      VEC3 phong = lights[l].color * pow(max(float(0), float(highlight.dot(-rayDir))), 10);

      shadowFlag = !getInShadow(point, (lights[l].position - point).normalized());
      c += shadowFlag * (sphereColors[hitID].cwiseProduct(max(float(0), l_norm) * lights[l].color + phong));    
    }
    pixelColor = c;
  } 
  
  if (cylinderHit) {
    c.setZero();
    for(int l = 0; l < lights.size(); l++) {
      VEC3 to_light = ((eye - point) + (lights[l].position - point)).normalized();

      float l_norm = normal.dot((lights[l].position - point).normalized());
      float h_norm = normal.dot(to_light);

      VEC3 highlight = ((-1 * (lights[l].position - point).normalized()) + 2 * l_norm * normal).normalized();
      VEC3 phong = lights[l].color * pow(max(float(0), float(highlight.dot(-rayDir))), 10);

      shadowFlag = !getInShadow(point, (lights[l].position - point).normalized());
      c += shadowFlag * boneColor.cwiseProduct(max(float(0), l_norm) * lights[l].color + phong);    
    }
    pixelColor = c;
  }
}
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
float clamp(float value)
{
  if (value < 0.0)      return 0.0;
  else if (value > 1.0) return 1.0;
  return value;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void renderImage(int& xRes, int& yRes, const string& filename) 
{
  // allocate the final image
  const int totalCells = xRes * yRes;
  float* ppmOut = new float[3 * totalCells];

  // compute image plane
  const float halfY = (lookingAt - eye).norm() * tan(45.0f / 360.0f * M_PI);
  const float halfX = halfY * 4.0f / 3.0f;

  const VEC3 cameraZ = (lookingAt - eye).normalized();
  const VEC3 cameraX = up.cross(cameraZ).normalized();
  const VEC3 cameraY = cameraZ.cross(cameraX).normalized();

  for (int y = 0; y < yRes; y++) 
    for (int x = 0; x < xRes; x++) 
    {
      const float ratioX = 1.0f - x / float(xRes) * 2.0f;
      const float ratioY = 1.0f - y / float(yRes) * 2.0f;
      const VEC3 rayHitImage = lookingAt + 
                               ratioX * halfX * cameraX +
                               ratioY * halfY * cameraY;
      const VEC3 rayDir = (rayHitImage - eye).normalized();

      // get the color
      VEC3 color;
      rayColor(eye, rayDir, color);

      // set, in final image
      ppmOut[3 * (y * xRes + x)] = clamp(color[0]) * 255.0f;
      ppmOut[3 * (y * xRes + x) + 1] = clamp(color[1]) * 255.0f;
      ppmOut[3 * (y * xRes + x) + 2] = clamp(color[2]) * 255.0f;
    }
  writePPM(filename, xRes, yRes, ppmOut);

  delete[] ppmOut;
}

//////////////////////////////////////////////////////////////////////////////////
// Load up a new motion captured frame
//////////////////////////////////////////////////////////////////////////////////
void setSkeletonsToSpecifiedFrame(int frameIndex)
{
  if (frameIndex < 0)
  {
    printf("Error in SetSkeletonsToSpecifiedFrame: frameIndex %d is illegal.\n", frameIndex);
    exit(0);
  }
  if (displayer.GetSkeletonMotion(0) != NULL)
  {
    int postureID;
    if (frameIndex >= displayer.GetSkeletonMotion(0)->GetNumFrames())
    {
      cout << " We hit the last frame! You might want to pick a different sequence. " << endl;
      postureID = displayer.GetSkeletonMotion(0)->GetNumFrames() - 1;
    }
    else 
      postureID = frameIndex;
    displayer.GetSkeleton(0)->setPosture(* (displayer.GetSkeletonMotion(0)->GetPosture(postureID)));
  }
}

void construct() { //construct static scene objects
  sphereCenters.push_back(VEC3(0.4,-0.2,1.4));
  sphereRadii.push_back(1.05);
  sphereColors.push_back(VEC3(0,0,1));

  sphereCenters.push_back(VEC3(0,-1500,0));
  sphereRadii.push_back(1500);
  sphereColors.push_back(VEC3(0.5,0.5,0.5));

  Triangle t1 = {
    VEC3(4, -0.1, 2), //a
    VEC3(2, -0.1, -2), //b
    VEC3(1, 2.5, -0.5), //c
    true
  };

  triangles.push_back(t1);
}

//////////////////////////////////////////////////////////////////////////////////
// Build a list of spheres in the scene
//////////////////////////////////////////////////////////////////////////////////
void buildScene()
{
  sphereCenters.clear();
  sphereRadii.clear();
  sphereColors.clear();
  cylinders.clear();
  displayer.ComputeBonePositions(DisplaySkeleton::BONES_AND_LOCAL_FRAMES);

  construct();

  // retrieve all the bones of the skeleton
  vector<MATRIX4>& rotations = displayer.rotations();
  vector<MATRIX4>& scalings  = displayer.scalings();
  vector<VEC4>& translations = displayer.translations();
  vector<float>& lengths     = displayer.lengths();

  // build a sphere list, but skip the first bone, 
  // it's just the origin
  int totalBones = rotations.size();
  for (int x = 1; x < totalBones; x++)
  {
    MATRIX4& rotation = rotations[x];
    MATRIX4& scaling = scalings[x];
    VEC4& translation = translations[x];

    VEC4 pelvisTranslation = translations[1];

    eye[0] = pelvisTranslation[0] - 4; 
    eye[1] = pelvisTranslation[1] + 0.5;
    eye[2] = pelvisTranslation[2];

    // get the endpoints of the cylinder
    VEC4 leftVertex(0,0,0,1);
    VEC4 rightVertex(0,0,lengths[x],1);

    leftVertex = rotation * scaling * leftVertex + translation;
    rightVertex = rotation * scaling * rightVertex + translation;

    // get the direction vector
    VEC3 direction = (rightVertex - leftVertex).head<3>();
    const float magnitude = direction.norm();
    direction *= 1.0 / magnitude;

    // how many spheres?
    const float sphereRadius = 0.05;
    const int totalSpheres = magnitude / (2.0 * sphereRadius);
    const float rayIncrement = magnitude / (float)totalSpheres;

    // store the spheres
    sphereCenters.push_back(leftVertex.head<3>());
    sphereRadii.push_back(0.05);
    sphereColors.push_back(VEC3(1,0,0));
    
    sphereCenters.push_back(rightVertex.head<3>());
    sphereRadii.push_back(0.05);
    sphereColors.push_back(VEC3(1,0,0));

    for (int y = 0; y < totalSpheres; y++)
    {
      VEC3 center = ((float)y + 0.5) * rayIncrement * direction + leftVertex.head<3>();

      //store the cylinders
      Cylinder tempCylinder = {
        (rightVertex + leftVertex).head<3>() / 2, //center point of cylinder??
        direction, //axis unit vector
        leftVertex.head<3>(),
        rightVertex.head<3>(),
        ((rightVertex - leftVertex).head<3>()).norm(), // height
        0.03,
      };

      cylinders.push_back(tempCylinder);

      // sphereCenters.push_back(center);
      // sphereRadii.push_back(0.05);
      // sphereColors.push_back(VEC3(1,0,0));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  string skeletonFilename("88.asf");
  string motionFilename("88_02.amc");
  // string skeletonFilename("02.asf");
  // string motionFilename("02_05.amc");
  
  // load up skeleton stuff
  skeleton = new Skeleton(skeletonFilename.c_str(), MOCAP_SCALE);
  skeleton->setBasePosture();
  displayer.LoadSkeleton(skeleton);

  // load up the motion
  motion = new Motion(motionFilename.c_str(), MOCAP_SCALE, skeleton);
  displayer.LoadMotion(motion);
  skeleton->setPosture(*(displayer.GetSkeletonMotion(0)->GetPosture(0)));

  // Note we're going 8 frames at a time, otherwise the animation
  // is really slow.
  Light light1 = {
    VEC3(0, 10, 0), 
    VEC3(0.5, 0.5, 0.5)
  };

  Light light2 = {
    VEC3(5, 10, -10), 
    VEC3(0.5, 0.5, 0.5)
  };

  Light light3 = {
    VEC3(-5, 10, 10), 
    VEC3(0.5, 0.5, 0.5)
  };

  lights.push_back(light1);
  lights.push_back(light2);
  lights.push_back(light3);

  for (int x = 0; x < 1200; x += 3)
  {
    setSkeletonsToSpecifiedFrame(x);
    buildScene();

    char buffer[256];
    snprintf(buffer, 256, "./frames/frame.%04i.ppm", x / 3);
    renderImage(windowWidth, windowHeight, buffer);
    cout << "Rendered " + to_string(x / 3) + " frames" << endl;
  }

  return 0;
}
