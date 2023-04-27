/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: irocha
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cmath>
#include <algorithm>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

#define PI 3.141529

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes (change back to 640 x 480)
#define WIDTH 320
#define HEIGHT 240

//the field of view of the camera
#define fov 60.0

double aspectRatio = (double) WIDTH / (double) HEIGHT;
double leftOfImage = -aspectRatio * tan((fov / 2) * PI / 180.0);
double rightOfImage = aspectRatio * tan((fov / 2) * PI / 180.0);
double bottomOfImage = -tan((fov / 2) * PI / 180.0);
double topOfImage = tan((fov / 2) * PI / 180.0);



unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

struct Ray {
    glm::vec3 originOfRay;
    glm::vec3 directionOfRay;
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

double triangleIntersect(Ray ray) {
    //Currently set to the max value of a double so as not to risk an actual value being above it and not replaced properly
    double triangleIntersectionValue = DBL_MAX;
    
    for (int i = 0; i < num_triangles; i++) {
        Triangle currTriangle = triangles[i];
        
        glm::vec3 A = glm::vec3(currTriangle.v[0].position[0], currTriangle.v[0].position[1], currTriangle.v[0].position[2]);
        glm::vec3 B = glm::vec3(currTriangle.v[1].position[0], currTriangle.v[1].position[1], currTriangle.v[1].position[2]);
        glm::vec3 C = glm::vec3(currTriangle.v[2].position[0], currTriangle.v[2].position[1], currTriangle.v[2].position[2]);
        glm::vec3 normal = glm::cross(B - A, C - A);

        bool isParallel = glm::dot(ray.directionOfRay, normal) == 0.0;
        
        //if the dot product is 0, there is no intersection between the ray and plane 
        if (!isParallel) {
            //From lecture 16 ray-polygon intersection II
            double t = (glm::dot(A, normal) - glm::dot(ray.originOfRay, normal)) / glm::dot(ray.directionOfRay, normal);
            //Check if point p is inside the triangle
            //t must be greater than 0 or else the intersection is behind the origin of the ray
            if (t > 0 && t < triangleIntersectionValue) {
                glm::vec3 P = ray.originOfRay + glm::vec3(ray.directionOfRay.x * t, ray.directionOfRay.y * t, ray.directionOfRay.z * t);
                //If the dot product is positive, it means that the point P is on the same side of the triangle as its normal vector, indicating that it is outside the triangle. 
                //If the dot product is negative, it means that the point P is on the opposite side of the triangle as its normal vector, indicating that it is inside the triangle.
                glm::vec3 AxP0 = glm::cross(B - A, P - A);
                glm::vec3 BxP0 = glm::cross(C - B, P - B);
                glm::vec3 CxP0 = glm::cross(A - C, P - C);
                if (glm::dot(AxP0, BxP0) > 0 && glm::dot(AxP0, CxP0) > 0) triangleIntersectionValue = t;
              
            }
        }

    }
    if (triangleIntersectionValue == DBL_MAX) return -1.0;
    return triangleIntersectionValue;
}

double sphereIntersect(Ray ray) {
    //Currently set to the max value of a double so as not to risk an actual value being above it and not replaced properly
    double sphereIntersectionValue = DBL_MAX;

    for (int i = 0; i < num_spheres; i++) {
        Sphere currSphere = spheres[i];
        //From lecture 16, Ray-Sphere Intersection II
        double b = 2.0 * (ray.directionOfRay.x * (ray.originOfRay.x - currSphere.position[0]) +
            ray.directionOfRay.y * (ray.originOfRay.y - currSphere.position[1]) +
            ray.directionOfRay.z * (ray.originOfRay.z - currSphere.position[2]));
        double c = pow(ray.originOfRay.x - currSphere.position[0], 2) +
            pow(ray.originOfRay.y - currSphere.position[1], 2) +
            pow(ray.originOfRay.z - currSphere.position[2], 2) -
            pow(currSphere.radius, 2);

        bool isRadicalPositive = pow(b, 2) - (4 * c) >= 0;
        if (isRadicalPositive) {
            double t0 = (-b + sqrt(pow(b, 2) - (4 * c))) / 2.0;
            double t1 = (-b - sqrt(pow(b, 2) - (4 * c))) / 2.0;
            if (t0 > 0 && t1 > 0) sphereIntersectionValue = std::min(t0, t1);
            else if (t0 > 0) sphereIntersectionValue = t0;
            else if (t1 > 0) sphereIntersectionValue = t1;
        }
    }
    if (sphereIntersectionValue == DBL_MAX) return -1.0;
    return sphereIntersectionValue;
}

void getColor(glm::vec3 position, glm::vec3 currColor) {

}
//MODIFY THIS FUNCTION
void draw_scene()
{
    float zeroArray[3] = { 0.0, 0.0, 0.0 };
    glm::vec3 camPos = glm::make_vec3(zeroArray);
    //currI and currJ used to move along image and draw each line
    double currI = leftOfImage;
    for (int i = 0; i < WIDTH; i++) {
        glPointSize(2.0);
        glBegin(GL_POINTS);

        double currJ = bottomOfImage;
        for (int j = 0; j < HEIGHT; j++) {
            //Will be updated with color values r, g, b
            glm::vec3 rayTracedColor = glm::vec3(0.0, 0.0, 0.0);
            double finalIntersectionValue;

            double currDirectionArray[3] = { currI, currJ, -1 };
            glm::vec3 currDirection = glm::make_vec3(currDirectionArray);
            glm::vec3 normalizedDirection = glm::normalize(currDirection);

            Ray emittedRay;
            emittedRay.originOfRay = camPos;
            emittedRay.directionOfRay = normalizedDirection;

            double triangleIntersectionValue = triangleIntersect(emittedRay);
            double sphereIntersectionValue = sphereIntersect(emittedRay);

            if (triangleIntersectionValue < 0 && sphereIntersectionValue < 0) {
                //There is no intersection, make white
                plot_pixel(i, j, 255.0, 255.0, 255.0);
            }
            else {
                //If it does not intersect with a triangle but intersects with a sphere, or it intersects with both but the sphere takes precedent
                if ((triangleIntersectionValue <= 0 && sphereIntersectionValue > 0) || (sphereIntersectionValue > 0 && sphereIntersectionValue < triangleIntersectionValue)) {
                    finalIntersectionValue = sphereIntersectionValue;
                    //temp make red to test intersection
                    rayTracedColor.x = 255.0;
                }
                //If it does not intersect with a sphere but intersects with a triangle, or it intersects with both but the triangle takes precedent

                else if ((sphereIntersectionValue <= 0 && triangleIntersectionValue > 0) || (triangleIntersectionValue > 0 && triangleIntersectionValue < sphereIntersectionValue)) {
                    finalIntersectionValue = triangleIntersectionValue;
                    //temp make green to test intersection
                    rayTracedColor.y = 255.0;
                }
                if (finalIntersectionValue > 0) plot_pixel(i, j, rayTracedColor.x, rayTracedColor.y, rayTracedColor.z);
            }

            currJ += (topOfImage - bottomOfImage) / (HEIGHT - 1);
        }
        glEnd();
        glFlush();
        currI += (rightOfImage - leftOfImage) / (WIDTH - 1);

    }
  
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

