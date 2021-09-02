#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include "SubLiME.h"

/*      B-spline data: knots, coefficients, degree          */

/* number of control point & degree of the spline */
#define CPN   8
#define DEG   3

/* knots of a uniform bspline */
float knots[CPN+DEG+1] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0 , 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};

/* number of enclosure segment over each control polygon segment */
#define SEG    1

/* slefe result */
#define SLFN ((CPN-3)*SEG+1)  // number of break points on the slefe
REAL upperx[SLFN], uppery[SLFN];
REAL lowerx[SLFN], lowery[SLFN];

/* control points */
#define sx  50
#define dx  60
#define sy  200
#define dy  60
REAL cx[CPN] = {sx, sx+dx, sx+2*dx, sx+3*dx, sx+4*dx, sx+5*dx, sx+6*dx, sx+7*dx};
REAL cy[CPN] = {sy, sy+dy, sy-2*dx, sy+5*dy, sy+5*dy, sy-2*dy, sy+dy, sy};

/* the point dragged by mouse */
int  MovingPoint =-1;
int  WindowWidth= 600, WindowHeight = 600;  // the window size

// openGL bspline handler 
GLUnurbsObj * line;
void draw_box(GLdouble px_i, GLdouble py_i, GLdouble mx_i, GLdouble my_i, 
          GLdouble px_i1, GLdouble py_i1, GLdouble mx_i1, GLdouble my_i1) ;

/* bound the input bspline */
void bound_bspline()
{
    // compute the upper and lower bounds for x and y components
    bspSlefe(cx, CPN, DEG, 1, SEG, upperx, lowerx, 1);
    bspSlefe(cy, CPN, DEG, 1, SEG, uppery, lowery, 1);
	// bspSlefe is a function provided by SubLiME library.
	//  refer to SubLiME.h for the description
}



void windowinit(void)
{
    /* set up viewing */
    /* 500 x 500 window with origin lower left */

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, WindowWidth, 0.0, WindowHeight, -WindowHeight, WindowHeight);

    glMatrixMode(GL_MODELVIEW); 

	InitBounds(); // initialization of SubLiME library *********** NEEDED

    bound_bspline();  // compute the bound for the input bspline

	line = gluNewNurbsRenderer();
	gluNurbsProperty(line, GLU_SAMPLING_TOLERANCE, 10);

}

 
void display( void )
{
	float ctrlpoints[CPN*3];

    /* temporary variable */
    int i;

    glClearColor(0.0, 0.0, 0.0, 0.0); /* white background */
    glClear(GL_COLOR_BUFFER_BIT);  /* clear the window */

    /* plots the control polygon */
    glColor3f(0.0, 1.0, 0.0);      /* light green */
    glLineWidth(3.0);       
    glBegin(GL_LINE_STRIP);
    for( i=0; i< CPN; i++)
    {
        glVertex2f(cx[i], cy[i]);
    }
    glEnd();
    glLineWidth(1.0);       

    /* plots the control polygon */

    glColor3f(0.0, 1.0, 0.0);      /* green */
    for( i=0; i< CPN; i++) {
		glPushMatrix();
		glTranslated(cx[i], cy[i], 1.0);
		glutSolidCube(8.0);
		glPopMatrix();
	}

    glColor3f(0.5, 0.5, 0.5);
    for(i=0;i<SLFN-1;i++) {
            int i1 = i+1;
            draw_box( upperx[i],  uppery[i],
                      lowerx[i],  lowery[i],
                      upperx[i1], uppery[i1],
                      lowerx[i1], lowery[i1]);
    }

    /* computes and plots the spline */
    glColor3f(1.0, 1.0, 1.0);      /* draw spline in white */
	for(i=0; i<CPN; i++)
	{
		ctrlpoints[i*3] = cx[i];
		ctrlpoints[i*3+1] = cy[i];
		ctrlpoints[i*3+2] = 0;
	}

	gluBeginCurve(line);
		gluNurbsCurve(line, DEG+CPN+1, knots, 3, ctrlpoints, (DEG+1),  GL_MAP1_VERTEX_3);
	gluEndCurve(line);

    glutSwapBuffers(); /* clear buffers */
}


/* mouse clicking */
void mouse(int button, int state, int x, int y)
{
    int i;

    x = (int)((float)x/WindowWidth*600);
    y = (int)((float)(WindowHeight-y)/WindowHeight*600);

    if(state==GLUT_DOWN)
    {
        /* find the point dragged by mouse */
        for(i=0;i<CPN;i++)
            if( (cx[i]-x)>-5.0 && (cx[i]-x)<5.0 &&
                (cy[i]-y)>-5.0 && (cy[i]-y)<5.0 )
            {
                MovingPoint = i;
            }
    }
    else if(state==GLUT_UP)
    {
        /* clear the point number */
        MovingPoint = -1;
    }
}

/* move the control point */
void motion(int x, int y)
{
    x = (int)((float)x/WindowWidth*600);
    y = (int)((float)(WindowHeight-y)/WindowHeight*600);

    if(MovingPoint >=0)  // moved control point
    {
        cx[MovingPoint] = x;
        cy[MovingPoint] = y;
		bound_bspline();
	    glutPostRedisplay();
    }
}

void keyboard(unsigned char key, int x, int y)
{
    if (key ==27)  /* press ESC to quit */
        exit(0);
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
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB); /* default, not needed */
    glutInitWindowSize(500,500); /* 500 x 500 pixel window */
    glutInitWindowPosition(0,0); /* place window top left on display */
    glutCreateWindow("SubLiME Package -- Bspline Curve Demo (SurfLab, UFL)"); /* window title */
    glutDisplayFunc(display); /* display callback invoked when window opened */
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(keyboard);
    glutReshapeFunc(reshape);
    windowinit(); /* set attributes */
    glutMainLoop(); /* enter event loop */
    return 0;
} 


/* draw a gray box indicating the enclosure */
void draw_box(GLdouble px_i, GLdouble py_i, GLdouble mx_i, GLdouble my_i, 
          GLdouble px_i1, GLdouble py_i1, GLdouble mx_i1, GLdouble my_i1) 
{
     glColor3f(0.5, 0.5, 0.5);
     glBegin(GL_POLYGON);
         glVertex2d(px_i,  py_i);
         glVertex2d(px_i1, py_i1);
         glVertex2d(px_i1, my_i1);
         glVertex2d(px_i,  my_i);
     glEnd();
     glBegin(GL_POLYGON);
         glVertex2d(mx_i,  py_i);
         glVertex2d(mx_i1, py_i1);
         glVertex2d(px_i1, py_i1);
         glVertex2d(px_i,  py_i);
     glEnd();
     glBegin(GL_POLYGON);
         glVertex2d(px_i,  my_i);
         glVertex2d(px_i1, my_i1);
         glVertex2d(mx_i1, my_i1);
         glVertex2d(mx_i,  my_i);
     glEnd();
     glBegin(GL_POLYGON);
         glVertex2d(mx_i,  my_i);
         glVertex2d(mx_i1, my_i1);
         glVertex2d(mx_i1, py_i1);
         glVertex2d(mx_i,  py_i);
     glEnd();

     glColor3f(1.0, 0.3, 0.3);
     glBegin(GL_LINE_LOOP);
         glVertex2d(px_i,  py_i); glVertex2d(mx_i,  py_i);
         glVertex2d(mx_i,  my_i); glVertex2d(px_i,  my_i);
     glEnd();
     glBegin(GL_LINE_LOOP);
         glVertex2d(px_i1,  py_i1); glVertex2d(mx_i1,  py_i1);
         glVertex2d(mx_i1,  my_i1); glVertex2d(px_i1,  my_i1);
     glEnd();
}     
