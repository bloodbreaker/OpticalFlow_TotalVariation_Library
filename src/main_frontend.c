#ifndef FRONTEND
#define FRONTEND

/*---------------------------------------------------------------------------*/
/*             copyright 2011 by Oliver Demetz and Andres Bruhn              */      
/*                               Mathematical Image Analysis Group           */
/*                               Saarland University                         */
/*                                                                           */
/* Please report errors to demetz@mia.uni-saarland.de                        */
/*                                                                           */
/* This C-Code will open a GLUT-window and display the image that has been   */
/* specified as commandline argument. By pressing the '.'-button the         */
/* computation is started (see compute()).                                   */
/*---------------------------------------------------------------------------*/

#include <GL/glut.h>
#include <math.h>

#include "frontend_lib.c"
#include "malloc_lib.c"
#include "io_lib.c"
#include "noise_lib.c"
#include "matrix_lib.c"
#include "arg_utils.c"

typedef unsigned char uchar;
                      
char   filename[80];    
float  ***u;                    // clean image
float  ***un;                   // noisy version of u
float  ***v;                    // filtering result      
uchar  ***p6;                   // uploaded to gpu
int    numchannels;             // 1 if pgm, 3 if ppm

int    nx,ny;                   // image dimensions
int    bx,by;                   // border sizes
float  hx,hy;                   // grid size

int    lastclickx,lastclicky;   // position of last mouse click
uchar  active_param;            // which option is active
int    num_images_side_by_side; // number of images shown side by side
GLuint gl_rgb_tex;              // texture identifyier  

float  sigma_noise;             // std. deviation of additive noise

/*---------------------------------------------------------------------------*/

void handleDraw()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

	glBindTexture(GL_TEXTURE_2D, gl_rgb_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, ny, num_images_side_by_side*nx, 0, 
				 GL_RGB, GL_UNSIGNED_BYTE, p6[0][0]);
	glBegin(GL_QUADS);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	glTexCoord2f(0, 0); glVertex3f(0,0,0);
	glTexCoord2f(0, 1); glVertex3f(num_images_side_by_side*nx,0,0);     
	glTexCoord2f(1, 1); glVertex3f(num_images_side_by_side*nx,ny,0);    
	glTexCoord2f(1, 0); glVertex3f(0,ny,0);       
	glEnd();
	glutSwapBuffers();
}

/*---------------------------------------------------------------------------*/

void writeParamsToStr(char* tgt)
{
	printTitleLine(tgt,"frontend  skeleton");
	showParamLinef(tgt,active_param,'N',"sigma noise",&sigma_noise);
	printFooter(tgt);
}

/*---------------------------------------------------------------------------*/

void showParams()
{
	char *str;
	str = (char*) calloc(2048,sizeof(char));
	writeParamsToStr(str);
	fputs(str,stdout);
	fflush(stdout);
	free(str);
}

/*---------------------------------------------------------------------------*/

void compose_image()
{    
	int   i,j,k;
	unsigned char val;
	for(k=0;k<numchannels;k++)
	for(i=bx; i < nx+bx;i++)
	for(j=by; j < ny+by;j++)	
	{
		val=float_to_uchar(v[k][i][j]);
		p6[i-bx][j-by][k]=val;
		if(numchannels==1) // if pgm-mode
			p6[i-bx][j-by][1]=p6[i-bx][j-by][2]=val;
	}
}

/*---------------------------------------------------------------------------*/

void compute()
{
	char *str;
	int  k;
	
	for(k=0;k<numchannels;k++)
	{
		add_gauss_noise(u[k],un[k],0,sigma_noise,nx,ny,bx,by);
		//compute here
		copy_matrix_2d(un[k],v[k],nx,ny,bx,by);
		//user clicked here
		v[k][lastclickx+bx][lastclicky+by]=255.0f;
	}
	
	// write new computation result to rotating sequence of files
	str = (char*) calloc(2048,sizeof(char));
	writeParamsToStr(str);
	write_pgm_image_comment(get_lastname(),nx,ny,bx,by,v[0],1,"%s",str);
	free(str);
	
	// display
	compose_image();
	glutPostRedisplay();
}

/*---------------------------------------------------------------------------*/

void handleKeyboardspecial(int key, int x, int y)
{   
	change_valf(key,active_param,'n',&sigma_noise);
	showParams();    
}

/*---------------------------------------------------------------------------*/

void handleKeyboard(unsigned char key, int x, int y)
{
	int h;
	switch (key)
	{
	case 46:  compute(); break;
	case 27:  exit(0);
	case 'o': h = glutGet(GLUT_WINDOW_HEIGHT);
			  glutReshapeWindow(h*num_images_side_by_side*nx/ny,h); 
			  glutPostRedisplay();
			  break;
	case '+': glutReshapeWindow(glutGet(GLUT_WINDOW_WIDTH)*1.5f,
								glutGet(GLUT_WINDOW_HEIGHT)*1.5f);
			  glutPostRedisplay();
			  break;
	case '-': glutReshapeWindow(glutGet(GLUT_WINDOW_WIDTH)*0.6666f,
								glutGet(GLUT_WINDOW_HEIGHT)*0.6666f);
			  glutPostRedisplay();
			  break;
	default:  active_param=key;
	}
	showParams();
}

/*---------------------------------------------------------------------------*/

void handleMouse(int button, int state, int cx, int cy)
{
	if(state == GLUT_DOWN)
	{
		//angles
		lastclickx = ((int)roundf((float)nx*num_images_side_by_side*cx
								/glutGet(GLUT_WINDOW_WIDTH)))%nx;
		lastclicky = (int)roundf((float)ny*cy/glutGet(GLUT_WINDOW_HEIGHT));
		glutPostRedisplay();
	}
}

/*---------------------------------------------------------------------------*/

void ReSizeGLScene(int Width, int Height)
{
	glViewport(0,0,Width,Height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho (0, num_images_side_by_side*nx, ny, 0, -1.0f, 1.0f);
	glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

/*---------------------------------------------------------------------------*/

int main (int argc, char* argv[])
{
	int maxgv;
	char * console_string = (char*)calloc(3000,sizeof(char));
	lastclickx=lastclicky=-1;
	active_param='\0';
	num_images_side_by_side=1;
	parse_arg_int   (argc,argv,"bordersizex",&bx,2,console_string);
	parse_arg_int   (argc,argv,"bordersizey",&by,2,console_string);
	parse_arg_float (argc,argv,"gridsizex",&hx,1.0f,console_string);
	parse_arg_float (argc,argv,"gridsizey",&hy,1.0f,console_string);
	parse_arg_string(argc,argv,"filename",filename,"test.ppm",console_string);
	parse_arg_float (argc,argv,"sigma_noise",&sigma_noise,150.0f,console_string);
	printf("------\n%s---------\n",console_string);
	if(filename[strlen(filename)-2]=='g')
	{
		numchannels=1;
		u=(float***)malloc(sizeof(float**));
		u[0] = read_pgm_image(filename,&nx,&ny,bx,by,&maxgv);
	}
	else if(filename[strlen(filename)-2]=='p')
	{
		numchannels=3;
		u = read_ppm_image(filename,&nx,&ny,bx,by,&maxgv);
	}

	calloc_multi(2,3,sizeof(float),numchannels,nx,ny,0,bx,by,0,0,0,&v,&un);
	calloc_multi(1,3,sizeof(unsigned char),num_images_side_by_side*nx,ny,
				 3,0,0,0, 0,0,0,&p6);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

	glutInitWindowSize(num_images_side_by_side*nx,ny);
	glutCreateWindow("Frontend Skeleton");

	glutDisplayFunc(handleDraw);
	glutKeyboardFunc(handleKeyboard);
	glutSpecialFunc(handleKeyboardspecial);
	glutMouseFunc(handleMouse);
	glutReshapeFunc(ReSizeGLScene);

	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &gl_rgb_tex);
	glBindTexture(GL_TEXTURE_2D, gl_rgb_tex);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	if(ny%4) glPixelStorei(GL_UNPACK_ALIGNMENT, 1) ;
	compute();
	showParams();
	glutMainLoop();

	return(0);
}


#endif

