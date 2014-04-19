#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <GL/glut.h>
#include <SOIL/SOIL.h>
#include <math.h>
#include <iostream>
#include <string.h>
#include <sstream>
#include "imageloader.h"
#include "cube.cpp"

using std::cerr;
using std::endl;
using std::string;

// Global variables
// Window/viewport
const int startwinsize = 400; // Starting window width & height, in pixels
const double neardist = 1.0;  // Near & far clipping distances
const double fardist = 10.0;
int flag1=0;

// Keyboard
const int ESCKEY = 27;        // ASCII value of escape character

// For scene
double angle1 = 0.0;           // Rotation amount
const double anglestep = 1.0; // Step for above
bool rotate = false;          // True if rotating
bool vertnormflag = true;     // true for vertex normals, false for facet normals
bool showwireframe = false;   // true if wire-frame outline is shown

// For Cube
const int CUBEDEPTH = 5;      // Depth for Cube; size is approx 2^depth
Cube<CUBEDEPTH> * thecube;    // The Cube itself; holds potentials, computes surf
Potentialtype thresh = 50;    // Current threshold
const Potentialtype threshchange = 2.;
                              // Amount to change above

// For metaballs
int flag2=0;
const double ballradius = 0.1;  // The radius at threshold 50
const double ballpositions[][3] = {  // x, y, z of centers. Should be in [-1,1]
   {  0,-.6, .6 },
   { .6,  0, .6 },
   { .6, .6,  0 },
   {  0, .6,-.6 },
   {-.6,  0,-.6 },
   {-.6,-.6,  0 },
   { 99, 99, 99 } };  // End mark

////////////////////////////////////////////////////////////////////////////////////////////////////
// handle_triangle
// Gets a triangle from marching cubes code.
// Draws it.
void handle_triangle(Triangle & t)
{
   glBegin(GL_TRIANGLES);
      if (!vertnormflag) glNormal3fv(t.facenorm);  // Normcoordtype is float
      for(int kk=0; kk<3; ++kk)
      {
         if (vertnormflag) glNormal3fv(t.norms[kk]);
         glVertex3fv(t.verts[kk]);  // Vertcoordtype is float
      }
   glEnd();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// tostring
// Convert argument to string class using operator<<
// Must include <sstream>
template<typename T>
std::string tostring(const T & input)
{
   std::ostringstream os;
   os << input;
   return os.str();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// printbitmap
// Prints the given string at the given raster position
//  using GLUT bitmap fonts.
// You probably don't want any rotations in the model/view
//  transformation when calling this function.
void printbitmap(const string msg, double x, double y)
{
   glRasterPos2d(x, y);
   for (string::const_iterator ii = msg.begin();
        ii != msg.end();
        ++ii)
   {
      glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *ii);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// idle
// The GLUT idle function
void idle()
{
   // Print OpenGL errors, if there are any (for debugging)
   if (GLenum err = glGetError())
   {
      cerr << "OpenGL ERROR: " << gluErrorString(err) << endl;
   }

   // Do rotation
   if (rotate)
   {
      angle1 += anglestep;
      glutPostRedisplay();
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// reshape
// The GLUT reshape function
void reshape(int w, int h)
{
   glViewport(0, 0, w, h);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluPerspective(50, double(w)/h, neardist, fardist);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// keyboard
// The GLUT keyboard function
void keyboard(unsigned char key, int x, int y)
{
   switch (key)
   {
        case ESCKEY:  // ESC: Quit
                exit(0);
                break;
        case 'r':     // R: Toggle rotation
        case 'R':
                rotate = !rotate;
                break;
        case 'n':     // N: Switch between vertex & facet normals
        case 'N':
                vertnormflag = !vertnormflag;
                glutPostRedisplay();
                break;
        case 'c':     // C: Toggle wire-frame outline
        case 'C':
                showwireframe = !showwireframe;
                glutPostRedisplay();
                break;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// special
// The GLUT special function
void special(int key, int x, int y)
{
   switch (key)
   {
   case GLUT_KEY_RIGHT:   // ->: Increase threshold
      thresh += threshchange;
      thecube->setthreshold(thresh);
      thecube->update();
      glutPostRedisplay();
      break;
   case GLUT_KEY_LEFT:  // <-: Decrease threshold
      thresh -= threshchange;
      thecube->setthreshold(thresh);
      thecube->update();
      glutPostRedisplay();
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// setupcube
// Sets up potentials in cube
// Potentials should range from 0 to 100 (?)
void setupcube()
{
   thecube->setthreshold(thresh);
   for (int ii=0; ii<thecube->POTSIDE; ++ii)
   {
      for (int jj=0; jj<thecube->POTSIDE; ++jj)
      {
         for (int kk=0; kk<thecube->POTSIDE; ++kk)
         {
            double x = double(ii)/(thecube->POTSIDE-1)*2.-1.;  // x, y, z in [-1,1]
            double y = double(jj)/(thecube->POTSIDE-1)*2.-1.;
            double z = double(kk)/(thecube->POTSIDE-1)*2.-1.;
            double f = 0.0;
            for (int ss=0; ballpositions[ss][0] < 90; ++ss)
            {
               // At its radius, ball adds 1/2 to potential (1/2 * 100 = 50).
               // Formula is a * (1/2a)^(distance/radius).
               const double a = 0.7;
               const double l2pla = log(2) + log(a);
               double dx = x-ballpositions[ss][0];
               double dy = y-ballpositions[ss][1];
               double dz = z-ballpositions[ss][2];
               double dist = sqrt(dx*dx+dy*dy+dz*dz);
               f += a * exp(dist * -l2pla/ballradius);
            }
            thecube->setpotential(ii, jj, kk, f*100);
         }
      }
   }
   thecube->update();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// init
// Initialization
// Called by main after window creation
void init()
{
   setupcube();

   // Misc OpenGL
   glEnable(GL_DEPTH_TEST);  // We're doing 3-D

   glClearColor(0.2, 0.2, 0.2, 1.0);  // Background color

   // Set up front material
   GLfloat diffuse_f[] = { 0.9, 0.2, 0.5, 1.0 };
   GLfloat specular_f[] = { 1.0, 1.0, 1.0, 1.0 };
   GLfloat shininess_f[] = { 100. };
   glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, diffuse_f);
   glMaterialfv(GL_FRONT, GL_SPECULAR, specular_f);
   glMaterialfv(GL_FRONT, GL_SHININESS, shininess_f);

   // Set up back material
   GLfloat diffuse_b[] = { 0.2, 0.7, 0.5, 1.0 };
   GLfloat specular_b[] = { 1.0, 1.0, 1.0, 1.0 };
   GLfloat shininess_b[] = { 100. };
   glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, diffuse_b);
   glMaterialfv(GL_BACK, GL_SPECULAR, specular_b);
   glMaterialfv(GL_BACK, GL_SHININESS, shininess_b);

   // Set up a light
   GLfloat light_position[] = { -1.0, 1.3, 1.0, 0.0 };
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);
   glEnable(GL_LIGHT0);

   // General lighting stuff
   glEnable(GL_LIGHTING);
   GLfloat lm_ambient[] = { 0.4, 0.4, 0.4, 1.0 };
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lm_ambient);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int flag=0;
enum MENU_TYPE
{
        MENU_ANGLE,
        MENU_SPEED,
        MENU_BALLS,
        MENU_QUIT,
        MENU_START,
};

////////////////////////////////////////////////////////////////////////////////////////////////////
const double PI = atan(1)*4;

GLUquadricObj *quadratic;
GLuint _textureId;          //The OpenGL id of the texture
int millisecs=20;
float maxangle=70;
float incrementmax = 6.5;
float angle = -maxangle;
bool clockwise = false;
int sphere=5;
int semueven = 1;
float diamm= 1.0;
float cubesradii = 0.125;

//Base                      ////size of wood base
float base_tamX = 7.5;
float base_tamY = 0.8;
float base_tamZ = 5.5;

float dist_esf_base = 1.0;  //Distance between the sphere and the base when the angle is 0
//Tubos
float tub_radio = 0.125;    //Radius of each cylinder
GLint slices = 32;          //Number of subdivisions around the z axis to draw each cylinder
GLint stacks = 32;          //Number of divisions along the z axis to draw each cylinder

float tub_tamX = 6.5;       // Size rectangular caramel
float tub_tamY = 5;      //It consists of pipes (including the diameter of each tube)
float tub_tamZ = 4.5;

// The length of the wire is calculated based on the size of the tubes,
// The diameter of the sphere and the distance between the sphere and the base
float largo_threads = tub_tamY-tub_radio-dist_esf_base-diamm-cubesradii/2;

//Camara
float camposx = 0.0;
float camposy = 0.5;
float camposz = 15.0;
float camrotx = 6.0;
float camroty = 40.0;
float camrotz = 0.0;

//axis
bool drawaxes = false;

////////////////////////////////////////////////////////////////////////////////////////////////////
GLuint loadTexture(Image* image) {
    GLuint textureId;
    glGenTextures(1, &textureId);
    glBindTexture(GL_TEXTURE_2D, textureId);
    glTexImage2D(GL_TEXTURE_2D,
                0,
                GL_RGB,
                image->width, image->height,
                0,
                GL_RGB,
                GL_UNSIGNED_BYTE,
                image->pixels);
    return textureId;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void initRendering() {
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glShadeModel(GL_SMOOTH);

    quadratic=gluNewQuadric();
    gluQuadricNormals(quadratic, GLU_SMOOTH);
    gluQuadricTexture(quadratic, GL_TRUE);

    Image* image = loadBMP("madera.bmp");
    _textureId = loadTexture(image);
    delete image;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void cameraposition() {
    glTranslatef(-camposx, -camposy, -camposz);
    glRotatef(camrotx, 1.0, 0.0, 0.0);
    glRotatef(camroty, 0.0, 1.0, 0.0);
    glRotatef(camrotz, 0.0, 0.0, 1.0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
double toDeg(double radian) {
    return radian*180/PI;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
double toRad(double degree) {
    return degree*PI/180;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int median(int n) {
     if(n%2==0)
        return n/2;
    else
        return (n+1)/2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void drawaxisss(void) {

    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);
        glVertex3f(-500.0,0.0,0.0);
        glVertex3f(500.0,0.0,0.0);
    glEnd();

    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_LINES);
        glVertex3f(0.0,-500.0,0.0);
        glVertex3f(0.0,500.0,0.0);
    glEnd();


    glColor3f(0.0, 0.0, 1.0);
    glBegin(GL_LINES);
        glVertex3f(0.0,0.0,-500.0);
        glVertex3f(0.0,0.0,500.0);
    glEnd();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void drawBase() {
    glPushMatrix();
    glTranslatef(0.0f, -(diamm/2+dist_esf_base+base_tamY/2), 0.0f);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, _textureId);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glColor3f(1.0f,1.0f,1.0f);

    glBegin(GL_QUADS);
        glNormal3f(0.0, 1.0f, 0.0f);
        glTexCoord2f(0.0f, 0.0f); glVertex3f(-base_tamX/2, base_tamY/2, base_tamZ/2); //Bottom Left of Texture & Plane
        glTexCoord2f(1.0f, 0.0f); glVertex3f(base_tamX/2, base_tamY/2, base_tamZ/2); // Bottom Right of Texture & Plane
        glTexCoord2f(1.0f, 1.0f); glVertex3f(base_tamX/2, base_tamY/2, -base_tamZ/2); // Top Right of Texture & Plane
        glTexCoord2f(0.0f, 1.0f); glVertex3f(-base_tamX/2, base_tamY/2, -base_tamZ/2); // Top Left of Texture & Plane
        glNormal3f(0.0, -1.0f, 0.0f);
        glTexCoord2f(0.0f, 0.0f); glVertex3f(-base_tamX/2, -base_tamY/2, base_tamZ/2);
        glTexCoord2f(1.0f, 0.0f); glVertex3f(base_tamX/2, -base_tamY/2, base_tamZ/2);
        glTexCoord2f(1.0f, 1.0f); glVertex3f(base_tamX/2, -base_tamY/2, -base_tamZ/2);
        glTexCoord2f(0.0f, 1.0f); glVertex3f(-base_tamX/2, -base_tamY/2, -base_tamZ/2);
        glNormal3f(-1.0, 0.0f, 0.0f);
        glTexCoord2f(0.0f, 0.0f); glVertex3f(-base_tamX/2, -base_tamY/2, -base_tamZ/2);
        glTexCoord2f(1.0f, 0.0f); glVertex3f(-base_tamX/2, -base_tamY/2, base_tamZ/2);
        glTexCoord2f(1.0f, 1.0f); glVertex3f(-base_tamX/2, base_tamY/2, base_tamZ/2);
        glTexCoord2f(0.0f, 1.0f); glVertex3f(-base_tamX/2, base_tamY/2, -base_tamZ/2);
        glNormal3f(-1.0, 0.0f, 0.0f);
        glTexCoord2f(0.0f, 0.0f); glVertex3f(base_tamX/2, -base_tamY/2, -base_tamZ/2);
        glTexCoord2f(1.0f, 0.0f); glVertex3f(base_tamX/2, -base_tamY/2, base_tamZ/2);
        glTexCoord2f(1.0f, 1.0f); glVertex3f(base_tamX/2, base_tamY/2, base_tamZ/2);
        glTexCoord2f(0.0f, 1.0f); glVertex3f(base_tamX/2, base_tamY/2, -base_tamZ/2);
        glNormal3f(0.0, 0.0f, 1.0f);
        glTexCoord2f(0.0f, 0.0f); glVertex3f(-base_tamX/2, -base_tamY/2, base_tamZ/2);
        glTexCoord2f(1.0f, 0.0f); glVertex3f(base_tamX/2, -base_tamY/2, base_tamZ/2);
        glTexCoord2f(1.0f, 1.0f); glVertex3f(base_tamX/2, base_tamY/2, base_tamZ/2);
        glTexCoord2f(0.0f, 1.0f); glVertex3f(-base_tamX/2, base_tamY/2, base_tamZ/2);
        glNormal3f(0.0, 0.0f, 1.0f);
        glTexCoord2f(0.0f, 0.0f); glVertex3f(-base_tamX/2, -base_tamY/2, -base_tamZ/2);
        glTexCoord2f(1.0f, 0.0f); glVertex3f(base_tamX/2, -base_tamY/2, -base_tamZ/2);
        glTexCoord2f(1.0f, 1.0f); glVertex3f(base_tamX/2, base_tamY/2, -base_tamZ/2);
        glTexCoord2f(0.0f, 1.0f); glVertex3f(-base_tamX/2, base_tamY/2, -base_tamZ/2);
    glEnd();

    glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void drawtuubes() {
    glPushMatrix();
    glTranslatef(0.0f, tub_tamY/2-diamm/2-dist_esf_base, 0.0f);
    glColor3f(0.69f, 0.69f, 0.69f);

    GLfloat mat_ambient[] = { 0.0, 0.0, 0.0, 0.0 };
    GLfloat mat_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 300.0 };
    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

    //Superior front
    glPushMatrix();
    glTranslatef(-tub_tamX/2, tub_tamY/2-tub_radio, tub_tamZ/2-tub_radio);
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    gluCylinder(quadratic,tub_radio,tub_radio,tub_tamX,slices, stacks);
    glPopMatrix();

    //Superior Posterior
    glPushMatrix();
    glTranslatef(-tub_tamX/2, tub_tamY/2-tub_radio, -(tub_tamZ/2-tub_radio));
    glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
    gluCylinder(quadratic,tub_radio,tub_radio,tub_tamX,slices, stacks);
    glPopMatrix();

    //left front
    glPushMatrix();
    glTranslatef(-(tub_tamX/2-tub_radio), tub_tamY/2-tub_radio, tub_tamZ/2-tub_radio);
    glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
    gluCylinder(quadratic,tub_radio,tub_radio,tub_tamY-tub_radio,slices, stacks);
    glPopMatrix();

    //left Posterior
    glPushMatrix();
    glTranslatef(-(tub_tamX/2-tub_radio), tub_tamY/2-tub_radio, -(tub_tamZ/2-tub_radio));
    glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
    gluCylinder(quadratic,tub_radio,tub_radio,tub_tamY-tub_radio,slices, stacks);
    glPopMatrix();

    //right front
    glPushMatrix();
    glTranslatef((tub_tamX/2-tub_radio), tub_tamY/2-tub_radio, tub_tamZ/2-tub_radio);
    glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
    gluCylinder(quadratic,tub_radio,tub_radio,tub_tamY-tub_radio,slices, stacks);
    glPopMatrix();

    //right Posterior
    glPushMatrix();
    glTranslatef((tub_tamX/2-tub_radio), tub_tamY/2-tub_radio, -(tub_tamZ/2-tub_radio));
    glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
    gluCylinder(quadratic,tub_radio,tub_radio,tub_tamY-tub_radio,slices, stacks);
    glPopMatrix();

    glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void tiedsphere(float angle) {
    glPushMatrix();

    glRotatef(angle, 0.0f, 0.0f, 1.0f);
    glTranslatef(0.0f, -largo_threads, 0.0f);

    glColor3f(0.8, 0.8, 0.8);


    GLfloat mat_ambient[] = { 0.0, 0.0, 0.0, 0.0 };
    GLfloat mat_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 300.0 };
    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);


    glutSolidCube(cubesradii);

    glPushMatrix();
    glTranslatef(0.0f, -(diamm/2), 0.0f);
    glutSolidSphere(diamm/2, 20, 20);
    glPopMatrix();

    glRotatef(-angle, 0.0f, 0.0f, 1.0f);

    //draw thread
    float distX = sin(toRad(angle))*largo_threads;
    float distY = cos(toRad(angle))*largo_threads;
    float distZ = tub_tamZ/2-tub_radio;

    glColor3f(0.72f, 0.54f, 0.0f);
    glBegin(GL_LINES);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(-distX, distY, -distZ);

        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(-distX, distY, distZ);
    glEnd();

    glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void setPoint(float x, float y){
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_POINTS);
        glVertex2f(x,y);
    glEnd();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void drawsphere() {
    glPushMatrix();
    glTranslatef(0.0f, tub_tamY-tub_radio-dist_esf_base-diamm/2, 0.0f);
    for(int i=1; i<=sphere; i++) {
        glPushMatrix();
        glTranslatef(-sphere/2.0f-diamm/2.0f+i*diamm,0.0f,0.0f);
        if(i<=semueven && angle<0)
            tiedsphere(angle);
        else if(i>sphere-semueven && angle>0)
            tiedsphere(angle);

        else
            tiedsphere(0);

        glPopMatrix();
    }

    glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void drawScene() {
    glClearColor(0.0f,0.0f,0.0f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    /** Iluminationn **/
    GLfloat ambientLight[] = {0.1f, 0.1f, 0.1f, 1.0f};
    GLfloat diffuseLight[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat specularLight[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat lightPos[] = {-camposx+tub_tamX, -camposy+tub_tamY, -camposz+tub_tamZ, 1.0f};
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);

    /** Viewing **/
    cameraposition();

    /** Modeling **/
    if(drawaxes)
        drawaxisss();
    drawBase();
    drawtuubes();
    drawsphere();
    glPointSize(20.0);
    //setPoint(a,b);
    glutSwapBuffers();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void update(int value) {
    float increment = incrementmax-abs(angle)/maxangle*incrementmax*0.85;

    if(clockwise and angle<=-maxangle) {
        clockwise=false;
    }
    else if(!clockwise and angle>=maxangle) {
        clockwise=true;
    }

    if(clockwise)
        angle -= increment;
    else
        angle += increment;

        printf("%d\n",maxangle);
        glutPostRedisplay();
        glutTimerFunc(millisecs, update, 0);

}

////////////////////////////////////////////////////////////////////////////////////////////////////
void handleSpecialKeys (int key, int x, int y) {
    int inc = 2.0;
    switch(key) {
            case GLUT_KEY_RIGHT:
            camroty += inc;
                glutPostRedisplay();
            break;

        case GLUT_KEY_LEFT:
            camroty -= inc;
            glutPostRedisplay();
            break;

        case GLUT_KEY_UP:
            camrotx -= inc;
            glutPostRedisplay();
            break;

        case GLUT_KEY_DOWN:
            camrotx += inc;
            glutPostRedisplay();
            break;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void handleResize(int w, int h) {
    /** Projection **/
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (float)w / (float)h, 1.0, 200.0);

    /** ViewPort **/
    glViewport(0, 0, w, h);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void handleMenuPrincipal(int m) {
    switch(m) {
        case 0:
            exit(0);
        case 1:
            drawaxes = !drawaxes;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void handleMenuTotalEsf(int m) {
    sphere = m-10;
    glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void handleMenuEnMovimiento(int m) {
    semueven = m-20;
    glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void menu(int item)
{
  printf( "in menu\n");
        switch (item)
        {
            case 'a':
                    sphere=3;
                    glutPostRedisplay();
                    break;
            case 'b':
                    sphere=5;
                    glutPostRedisplay();
                    break;
            case 'c':
                    sphere=7;
                    glutPostRedisplay();
                    break;

           case 'd':
                    if(angle<=45 || angle>=-45)
                    {
                        maxangle=angle=45;
                        glutPostRedisplay();
                    }
                    break;
            case 'e':
                    if(angle<=60 || angle>=-60)
                    {
                        maxangle=angle=60;
                        glutPostRedisplay();
                    }
                    break;
            case 'f':
                   if(angle<=90 || angle>=-90)
                    {
                        maxangle=angle=90;
                        glutPostRedisplay();
                    }
                    break;
            case 'g':
                    millisecs=50;
                    glutPostRedisplay();
                    break;
            case 'h':
                    millisecs=25;
                    glutPostRedisplay();
                    break;
            case 'i':
                    millisecs=10;
                    glutPostRedisplay();
                    break;

            case 's':
                    if('s')
                        glutDisplayFunc(drawScene);
                    else
                        glutPostRedisplay();
                    break;

            case 'q':
                    if('q')
                        exit(0);
                    else
                        glutPostRedisplay();
                    break;
        }
        return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void *currentfont;
GLuint tex_2d, tex_2d1, tex_2d2;

////////////////////////////////////////////////////////////////////////////////////////////////////
void drawstring( int x, int y, char *st) //this function will print the character
{
	int l,i;
	l=strlen( st ); // see how many characters are in text string.
	glRasterPos2i( x, y); // location to start printing text
	for( i=0; i < l; i++)  // loop until i is greater then l
	{

		glutBitmapCharacter(currentfont, st[i]);// Print a character on the screen
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void setfont(void *font)
{
    currentfont=font;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// display
// The GLUT display function
void display()
{
     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_LIGHTING);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
      glLoadIdentity();
      gluOrtho2D(-1., 1., -1., 1.);
      glColor3d(1.0, 1.0, 1.0);
      printbitmap("Metaballs Demo", -0.9, 0.9);
      printbitmap("R      Toggle rotation", -0.9, 0.8);
      printbitmap("N      Switch vertex <-> facet normals", -0.9, 0.7);
      printbitmap("C      Toggle wire-frame cube", -0.9, 0.6);
      printbitmap("<- ->  Change threshold [" + tostring(thresh) + "]", -0.9, 0.5);
      printbitmap("ESC    Quit", -0.9, 0.4);
   glPopMatrix();
   glMatrixMode(GL_MODELVIEW);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_LIGHTING);

      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glTranslated(0., 0., -4.);
      glRotated(angle1, 0.5, 1.0, 0.0);

      // Draw surface
      thecube->surfspit(handle_triangle);

      // Draw wire-frame outline, if shown
      if (showwireframe)
      {
         glDisable(GL_LIGHTING);
         glColor3d(1.0, 1.0, 0.0);
         glBegin(GL_LINE_LOOP);
            glVertex3d(-1.0, -1.0, -1.0);
            glVertex3d(1.0, -1.0, -1.0);
            glVertex3d(1.0, 1.0, -1.0);
            glVertex3d(-1.0, 1.0, -1.0);
         glEnd();
         glBegin(GL_LINE_LOOP);
            glVertex3d(-1.0, -1.0, 1.0);
            glVertex3d(1.0, -1.0, 1.0);
            glVertex3d(1.0, 1.0, 1.0);
            glVertex3d(-1.0, 1.0, 1.0);
         glEnd();
         glBegin(GL_LINES);
            glVertex3d(-1.0, -1.0, -1.0);
            glVertex3d(-1.0, -1.0, 1.0);

            glVertex3d(1.0, -1.0, -1.0);
            glVertex3d(1.0, -1.0, 1.0);

            glVertex3d(1.0, 1.0, -1.0);
            glVertex3d(1.0, 1.0, 1.0);

            glVertex3d(-1.0, 1.0, -1.0);
            glVertex3d(-1.0, 1.0, 1.0);
         glEnd();

      }
   glPopMatrix();
   glutSwapBuffers();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void nextp2()
{
    thecube = new Cube<CUBEDEPTH>(thresh, 0, 2, -1, -1, -1); assert(thecube != NULL);

   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

   // Make a window
   glutInitWindowSize(700,700);
   glutInitWindowPosition(50, 50);
   glutCreateWindow("Metaballs");

   // Initialize GL states & register callbacks
   init();
   glutDisplayFunc(display);
   glutIdleFunc(idle);
   glutReshapeFunc(reshape);
   glutKeyboardFunc(keyboard);
   glutSpecialFunc(special);

}

////////////////////////////////////////////////////////////////////////////////////////////////////
void nextp1()
{
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 600);

    glutCreateWindow("newtons cradle");
    initRendering();

    glutDisplayFunc(drawScene);

    int MENU_BALLS = glutCreateMenu(menu);
    glutAddMenuEntry("3", 'a');
    glutAddMenuEntry("5", 'b');
    glutAddMenuEntry("7", 'c');

    int MENU_ANGLE = glutCreateMenu(menu);
    glutAddMenuEntry("45 deg", 'd');
    glutAddMenuEntry("70 deg", 'e');
    glutAddMenuEntry("100 deg", 'f');

    int MENU_SPEED = glutCreateMenu(menu);
    glutAddMenuEntry("slow", 'g');
    glutAddMenuEntry("medium", 'h');
    glutAddMenuEntry("fast", 'i');

    int MENU_QUIT = glutCreateMenu(menu);
    glutAddMenuEntry("Yes", 'q');
    glutAddMenuEntry("No", 'n');

    glutCreateMenu(menu);

    glutAddSubMenu("No. of balls", MENU_BALLS);
    glutAddSubMenu("Angle of rotation", MENU_ANGLE);
    glutAddSubMenu("Movement speed", MENU_SPEED);
    //glutAddSubMenu("Wanna Start??", MENU_START);
    glutAddSubMenu("Quit??", MENU_QUIT);

    glutAttachMenu(GLUT_RIGHT_BUTTON);

    glutReshapeFunc(handleResize);
    glutSpecialFunc(handleSpecialKeys);
    glutTimerFunc(millisecs, update, 0);

}

////////////////////////////////////////////////////////////////////////////////////////////////////
void init1()
{
    glColor3f(0.0,0.0,0.0);
    glPointSize(5.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0,1200.0,0.0,700.0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void frontpg()
{
	glClear(GL_COLOR_BUFFER_BIT);
    glClearColor(1.0, 1.0, 1.0, 1.0);

    /*back*/
    glEnable(GL_TEXTURE_2D);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	tex_2d = SOIL_load_OGL_texture
		 (	 "bk.jpg",
			 SOIL_LOAD_AUTO,
			 SOIL_CREATE_NEW_ID,
			 SOIL_FLAG_COMPRESS_TO_DXT
		 );
	glBindTexture(GL_TEXTURE_2D, tex_2d);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glBegin(GL_POLYGON);
		glTexCoord2f(1.0, 0.0);
		glVertex2f(0,0);
		glTexCoord2f(0.0, 0.0);
		glVertex2f(1200,0);
		glTexCoord2f(0.0, 1.0);
		glVertex2f(1200,1200);
		glTexCoord2f(1.0, 1.0);
		glVertex2f(0,1200);
	glEnd();

    setfont(GLUT_BITMAP_TIMES_ROMAN_24);
    glColor3f(0.0,0.0,0.0);
    drawstring(200.0,647.0,"```````````````````````````````````````````````````````````````````````````````````````````````````````````````````");
    drawstring(200.0,530.0,"```````````````````````````````````````````````````````````````````````````````````````````````````````````````````");

    glColor3f(1.0,0.0,0.0);
    drawstring(450.0,600.0, "PESIT Bangalore South Campus");
    drawstring(370.0,575.0, "Department of Computer Science & Engineering");

    glColor3f(0.0,0.0,1.0);
    setfont(GLUT_BITMAP_HELVETICA_18);
    drawstring(200.0,123.0, "Submitted By,");
    drawstring(230.0,95.0, "-Jujare Vinayaka(1PE11CS039)");
    drawstring(230.0,67.0, "-Manikandan.R(1PE11CS046)");
    drawstring(900.0,120.0, "Guide,");
   	drawstring(930.0,90.0, "-Prof.Shubha Raj.K.B");

   	setfont(GLUT_BITMAP_TIMES_ROMAN_24);
 	glColor3f(0.0,1.0,0.0);
    drawstring(450.0,480.0, "A PROJECT WORK 0N");

    glColor3f(1.0,1.0,1.0);
    setfont(GLUT_BITMAP_TIMES_ROMAN_24);
    glColor3f(1.0,0.7,0.0);
    drawstring(230.0, 440.0, "** COMPUTER GRAPHICS AND VISUALIZATION LABORATORY **");
    drawstring(530,410," (10CSL67)");

    glColor3f(1.0,0.0,0.0);
    drawstring(220.0,360.0, "  TITLE    :     NEWTONS CRADLE      &      METABALLS SIMULATION");

    glColor3f(1.0,1.0,1.0);

    glColor3f(1.0,0.0,0.0);
    drawstring(420.0,320.0, "press 'n'");

    /*newtons cradle*/
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	tex_2d1 = SOIL_load_OGL_texture
		 (
			 "newton.jpg",
			 SOIL_LOAD_AUTO,
			 SOIL_CREATE_NEW_ID,
			 SOIL_FLAG_COMPRESS_TO_DXT
		 );
	glBindTexture(GL_TEXTURE_2D, tex_2d1);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glBegin(GL_POLYGON);
		glTexCoord2f(1.0, 1.0);
		glVertex2f(400,180);
		glTexCoord2f(0.0, 1.0);
		glVertex2f(520,180);
		glTexCoord2f(0.0, 0.0);
		glVertex2f(520,300);
		glTexCoord2f(1.0, 0.0);
		glVertex2f(400,300);
	glEnd();

    glColor3f(1.0,0.0,0.0);
    drawstring(760.0,320.0, "press 'm'");

	/*metaballs*/
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	tex_2d1 = SOIL_load_OGL_texture
		 (
			 "metaball.png",
			 SOIL_LOAD_AUTO,
			 SOIL_CREATE_NEW_ID,
			 SOIL_FLAG_COMPRESS_TO_DXT
		 );
	glBindTexture(GL_TEXTURE_2D, tex_2d1);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glBegin(GL_POLYGON);
		glTexCoord2f(1.0, 1.0);
		glVertex2f(740,180);
		glTexCoord2f(0.0, 1.0);
		glVertex2f(860,180);
		glTexCoord2f(0.0, 0.0);
		glVertex2f(860,300);
		glTexCoord2f(1.0, 0.0);
		glVertex2f(740,300);
	glEnd();

	/*pes logo*/
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	tex_2d1 = SOIL_load_OGL_texture
		 (
			 "logo1.png",
			 SOIL_LOAD_AUTO,
			 SOIL_CREATE_NEW_ID,
			 SOIL_FLAG_COMPRESS_TO_DXT
		 );
	glBindTexture(GL_TEXTURE_2D, tex_2d1);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glBegin(GL_POLYGON);
		glTexCoord2f(1.0, 1.0);
		glVertex2f(1160,550);
		glTexCoord2f(0.0, 1.0);
		glVertex2f(1085,550);
		glTexCoord2f(0.0, 0.0);
		glVertex2f(1085,650);
		glTexCoord2f(1.0, 0.0);
		glVertex2f(1160,650);
	glEnd();

	/*cse*/
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	tex_2d2 = SOIL_load_OGL_texture
		 (
			 "cse.png",
			 SOIL_LOAD_AUTO,
			 SOIL_CREATE_NEW_ID,
			 SOIL_FLAG_COMPRESS_TO_DXT
		 );
	glBindTexture(GL_TEXTURE_2D, tex_2d2);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glBegin(GL_POLYGON);
		glTexCoord2f(1.0, 1.0);
		glVertex2f(125,550);
		glTexCoord2f(0.0, 1.0);
		glVertex2f(50,550);
		glTexCoord2f(0.0, 0.0);
		glVertex2f(50,650);
		glTexCoord2f(1.0, 0.0);
		glVertex2f(125,650);
	glEnd();
    glDisable(GL_TEXTURE_2D);
    glutSwapBuffers();
    glFlush();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void mydisplay(void)
{
    glClear(GL_COLOR_BUFFER_BIT);
    if(flag==0)
        frontpg();
    if(flag==1)
        nextp1();
    if(flag1==1)
        nextp2();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void myKeyboardFunc( unsigned char key, int x, int y )
{
    switch(key)
    {
        case 'n':if(flag==0)
                flag=1;
                break;
        case 'm':if(flag1==0)
                flag1=1;
                break;
    }
    mydisplay();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1200,700);
    glutInitWindowPosition(0,0);
    glutCreateWindow("WELCOME");
    glutKeyboardFunc(myKeyboardFunc);
    glutDisplayFunc(frontpg);
    init1();
    glutMainLoop();
    return 0;
}
