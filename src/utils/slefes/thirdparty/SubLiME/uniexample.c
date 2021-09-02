#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include "SubLiME.h"

#define DEG 7            // degree
#define CPNUM (DEG + 1)  // number of control points

#define MAXSEG 9          // number of segments
#define PTS (MAXSEG + 1)  // number of breaking points

#define DIM 3  // dimension

#define REAL double  // data type

// input ---------------------------------
//
GLdouble ctrlpoints[CPNUM * DIM] = {
  30.0,  40.0,  0.0,  // control point 1
  130.0, 340.0, 0.0,  // control point 2
  230.0, 340.0, 0.0,  // control point 3
  330.0, 40.0,  0.0,  // control point 4
  430.0, 40.0,  0.0,  // control point 5
  530.0, 40.0,  0.0,  // control point 6
};

int seg = 5;

// results -------------------------------------------------
//
GLdouble Upper[PTS * DIM];
GLdouble Lower[PTS * DIM];

// function prototypes ------------------------------------
//

// draw a gray box
void draw_box(GLdouble px_i, GLdouble py_i, GLdouble mx_i, GLdouble my_i, GLdouble px_i1, GLdouble py_i1,
              GLdouble mx_i1, GLdouble my_i1);

// Some global variables -------------------------------------------
//
int MovingPoint = -1; /* the point dragged by mouse */
int WindowWidth = 600, WindowHeight = 600;

// code starts ---------------------------------------------------------
//

// this function demonstrates how to use the uniSlefe() function
//
// compute the upper and lower slefe for the input function
void bound_bezier()
{
  int m;

  // compute slefe for each component
  for (m = 0; m < DIM; m++)
    uniSlefe(&ctrlpoints[m], DIM, DEG, seg, &Upper[m], &Lower[m], DIM);
  // uniSlefe is a function provided by SubLiME library.
  //  refer to SubLiME.h for the description
}

// ---------------------------------------------------------------------

// Rest of the code enable the user to click on the control points
// and modify the curve and enclosures interactively

// glut window initialization
void windowinit(void)
{
  /* set up viewing */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, WindowWidth, 0.0, WindowHeight, -WindowHeight, WindowHeight);
  glMatrixMode(GL_MODELVIEW);

  // initialize the bounds
  InitBounds();  // initialization of SubLiME library *********** NEEDED

  bound_bezier();  // compute the slefe path for the input function
}

/* display function */
void display(void)
{
  int i, i1;

  // glClearColor(1.0, 1.0, 1.0, 0.0); /* white background */
  glClearColor(0.0, 0.0, 0.0, 0.0); /* white background */
  glClear(GL_COLOR_BUFFER_BIT);     /* clear the window */

  // plot the enclosures
  //
  // component enclosures are drawn with all the combinations
  glColor3f(0.5, 0.5, 0.5);
  for (i = 0; i < seg; i++)
  {
    i1 = i + 1;
    draw_box(Upper[i * DIM + 0], Upper[i * DIM + 1], Lower[i * DIM + 0], Lower[i * DIM + 1], Upper[i1 * DIM + 0],
             Upper[i1 * DIM + 1], Lower[i1 * DIM + 0], Lower[i1 * DIM + 1]);
  }

  /* Control polygon */
  glColor3f(0.0, 1.0, 0.0); /* green */
  glLineWidth(3.0);
  glBegin(GL_LINE_STRIP);
  for (i = 0; i < CPNUM; i++)
  {
    glVertex2f(ctrlpoints[i * DIM + 0], ctrlpoints[i * DIM + 1]);
  }
  glEnd();
  glLineWidth(1.0);

  /* Control points */
  glColor3f(0.0, 1.0, 0.0); /* green */
  for (i = 0; i < CPNUM; i++)
  {
    glPushMatrix();
    glTranslated(ctrlpoints[i * DIM + 0], ctrlpoints[i * DIM + 1], 1.0);
    glutSolidCube(8.0);
    glPopMatrix();
  }

  /* The Bezier curve */
  glColor3f(1.0, 1.0, 0.0); /* draw spline in Yellow */
  glLineWidth(1.0);
  glMap1d(GL_MAP1_VERTEX_3, 0.0, 1.0, DIM, DEG + 1, ctrlpoints);
  glEnable(GL_MAP1_VERTEX_3);
  glBegin(GL_LINE_STRIP);
  for (i = 0; i <= 50; i++)
    glEvalCoord1d((GLfloat)i / 50);
  glEnd();

  glutSwapBuffers(); /* clear buffers */
}

/* click a control point can move that point */
void mouse(int button, int state, int x, int y)
{
  int i;

  x = (int)((float)x / WindowWidth * 600);
  y = (int)((float)(WindowHeight - y) / WindowHeight * 600);

  if (state == GLUT_DOWN)
  {
    /* find the point dragged by mouse */
    for (i = 0; i < CPNUM; i++)
    {
      if ((ctrlpoints[i * DIM + 0] - x) > -5.0 && (ctrlpoints[i * DIM + 0] - x) < 5.0 &&
          (ctrlpoints[i * DIM + 1] - y) > -5.0 && (ctrlpoints[i * DIM + 1] - y) < 5.0)
      {
        MovingPoint = i;
      }
    }
  }
  else if (state == GLUT_UP)
  {
    /* clear the point number */
    MovingPoint = -1;
  }
}

/* mouse movement change the position of the picked point */
void motion(int x, int y)
{
  x = (int)((float)x / WindowWidth * 600);
  y = (int)((float)(WindowHeight - y) / WindowHeight * 600);

  if (MovingPoint >= 0)  // moved control point
  {
    ctrlpoints[MovingPoint * DIM + 0] = x;
    ctrlpoints[MovingPoint * DIM + 1] = y;

    bound_bezier();  // recompute the slefe
    glutPostRedisplay();
  }
}

void keyboard(unsigned char key, int x, int y)
{
  if (key == 27) /* press ESC to quit */
    exit(0);
  if (key == '+' || key == '=')
  {
    if (seg < MAXSEG)
      seg++;
    bound_bezier();
    glutPostRedisplay();
  }
  if (key == '_' || key == '-')
  {
    if (seg > 2)
      seg--;
    bound_bezier();
    glutPostRedisplay();
  }
}

void reshape(int width, int height)
{
  WindowWidth = width;
  WindowHeight = height;

  glViewport(0, 0, width, height);
  glutReshapeWindow(width, height);
  glutPostRedisplay();
}

int main(int argc, char** argv)
{
  /* Standard GLUT initialization */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);                             /* default, not needed */
  glutInitWindowSize(WindowWidth, WindowHeight);                           /* 600 x 600 pixel window */
  glutInitWindowPosition(0, 0);                                            /* place window top left on display */
  glutCreateWindow("SubLiME Package -- Bezier Curve Demo (SurfLab, UFL)"); /* window title */
  glutDisplayFunc(display); /* display callback invoked when window opened */
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutKeyboardFunc(keyboard);
  glutReshapeFunc(reshape);

  windowinit();   /* set attributes */
  glutMainLoop(); /* enter event loop */

  return 0;
}

/* draw a gray box indicating the enclosure */
void draw_box(GLdouble px_i, GLdouble py_i, GLdouble mx_i, GLdouble my_i, GLdouble px_i1, GLdouble py_i1,
              GLdouble mx_i1, GLdouble my_i1)
{
  glBegin(GL_POLYGON);
  glVertex2d(px_i, py_i);
  glVertex2d(px_i1, py_i1);
  glVertex2d(px_i1, my_i1);
  glVertex2d(px_i, my_i);
  glEnd();
  glBegin(GL_POLYGON);
  glVertex2d(mx_i, py_i);
  glVertex2d(mx_i1, py_i1);
  glVertex2d(px_i1, py_i1);
  glVertex2d(px_i, py_i);
  glEnd();
  glBegin(GL_POLYGON);
  glVertex2d(px_i, my_i);
  glVertex2d(px_i1, my_i1);
  glVertex2d(mx_i1, my_i1);
  glVertex2d(mx_i, my_i);
  glEnd();
  glBegin(GL_POLYGON);
  glVertex2d(mx_i, my_i);
  glVertex2d(mx_i1, my_i1);
  glVertex2d(mx_i1, py_i1);
  glVertex2d(mx_i, py_i);
  glEnd();
  glColor3f(1.0, 0.3, 0.3);
  glBegin(GL_LINE_LOOP);
  glVertex2d(px_i, py_i);
  glVertex2d(mx_i, py_i);
  glVertex2d(mx_i, my_i);
  glVertex2d(px_i, my_i);
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex2d(px_i1, py_i1);
  glVertex2d(mx_i1, py_i1);
  glVertex2d(mx_i1, my_i1);
  glVertex2d(px_i1, my_i1);
  glEnd();
  glColor3f(0.5, 0.5, 0.5);
}

//  code ends -----------------------------------------------
