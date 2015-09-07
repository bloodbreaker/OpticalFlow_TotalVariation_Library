#ifndef FRONTEND
#define FRONTEND

/*-------------------------------------------------------------------------------------*/
/*             copyright 2011 by Oliver Demetz and Andres Bruhn              */      
/*                               Mathematical Image Analysis Group           */
/*                               Saarland University                         */
/*                                                                           */
/* Please report errors to demetz@mia.uni-saarland.de                        */
/*                                                                           */
/* This C-Code will open a GLUT-window and display the image that has been   */
/* specified as commandline argument. By pressing the '.'-button the         */
/* computation is started (see compute()).                                   */
/*-------------------------------------------------------------------------------------*/

#include <GL/glut.h>
#include <math.h>

#include "frontend_lib.c"
#include "malloc_lib.c"
#include "io_lib.c"
#include "noise_lib.c"
#include "matrix_lib.c"
#include "arg_utils.c"
#include "color_lib.c"
#include "horn_schunck_warp.c"
#include "of_lib.c"
#include "diffusivity_lib.c"

typedef unsigned char uchar;
char basename[200];
float  **u;                     // optic flow field x-component
float  **v;                     // optic flow field y-component
float  **u_truth;               // ground truth
float  **v_truth;               // ground truth
float  ***of_rgb;                // color plot of flow field
float  **f1;                    // first  image
float  **f2;                    // second image
uchar  ***p6;                   // uploaded to gpu
int    nx,ny;                   // image dimensions
int    bx,by;                   // border sizes

float  alpha;                   // smoothness weight
float  epsilon_d;               // diffusivity param data term
int    diffusivity_type_d;      // nonlinearity function data term
int    num_iterations_inner;    // number of solver iterations
int    num_iterations_outer;    // number of fixed point iterations (nonlinearity)
float  omega;                   // SOR overrelaxation parameter
float  warp_eta;                // warping reduction factor between levels
int    max_warp_levels;         // maximal number of warping levels
float  max_displacement;        // max length of flowfield that is shown

float  aae;                     // Average Angular Error
float  al2e;                    // Average Endpoint Error (L2)
int    lastclickx,lastclicky;   // position of last mouse click
uchar  active_param;            // which option is active
int    num_images_side_by_side; // number of images shown side by side
GLuint gl_rgb_tex;              // texture identifyier  

/*-------------------------------------------------------------------------------------*/

void handleDraw()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

	glBindTexture(GL_TEXTURE_2D, gl_rgb_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, ny, num_images_side_by_side*nx, 0, 
				 GL_RGB, GL_UNSIGNED_BYTE, &(***p6));
	glBegin(GL_QUADS);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	glTexCoord2f(0, 0); glVertex3f(0,0,0);
	glTexCoord2f(0, 1); glVertex3f(num_images_side_by_side*nx,0,0);     
	glTexCoord2f(1, 1); glVertex3f(num_images_side_by_side*nx,ny,0);    
	glTexCoord2f(1, 0); glVertex3f(0,ny,0);       
	glEnd();
	glutSwapBuffers();
}

/*-------------------------------------------------------------------------------------*/

void writeParamsToStr(char* tgt)
{
	printTitleLine(tgt,"optic flow frontend");
	showParamLinef(tgt,active_param,'a',"Smoothness weight alpha",&alpha);
	showParamLinef(tgt,active_param,'E',"Diffusivity param Data",&epsilon_d);
	showParamLinei(tgt,active_param,'i',"Number inner iterations",&num_iterations_inner);
	showParamLinei(tgt,active_param,'I',"Number outer iterations",&num_iterations_outer);
	showParamLineStr(tgt,active_param,'d',"Diffusivity type Data",diffusivity_names[diffusivity_type_d]);
	showParamLinef(tgt,active_param,'O',"Overrelaxation omega",&omega);
	showParamLinef(tgt,active_param,'e',"Warping eta",&warp_eta);
	showParamLinei(tgt,active_param,'m',"Maximum number warping levels",&max_warp_levels);
	showParamLinef(tgt,active_param,'l',"Maximal length displacement",&max_displacement);
	printFooter(tgt);
	showParamLinef(tgt,active_param,' ',"Average Angular  Error",&aae);
	showParamLinef(tgt,active_param,' ',"Average Endpoint Error",&al2e);
	printFooter(tgt);
}

/*-------------------------------------------------------------------------------------*/

void showParams()
{
	char *str;
	str = (char*) calloc(2048,sizeof(char));
	writeParamsToStr(str);
	fputs(str,stdout);
	fflush(stdout);
	free(str);
}

/*-------------------------------------------------------------------------------------*/
// fill texture p6 for display
void compose_image()
{    
	int   i,j;
	flow_field_to_image(u,v,of_rgb,max_displacement,nx,ny,bx,by);
	for(i=bx; i < nx+bx;i++)
	for(j=by; j < ny+by;j++)	
	{
		// grey value image1
		uchar *tgt = p6[i-bx][j-by];
		*tgt=*(tgt+1)=*(tgt+2)=float_to_uchar(f1[i][j]);
		// flow field as color plot
		p6[i-bx+nx][j-by][0] = (uchar)of_rgb[0][i][j];
		p6[i-bx+nx][j-by][1] = (uchar)of_rgb[1][i][j];
		p6[i-bx+nx][j-by][2] = (uchar)of_rgb[2][i][j];
	}
}

/*-------------------------------------------------------------------------------------*/

void compute()
{
	char *str;
	
	/* compute from function	*/
	
	HORN_SCHUNCK_MAIN(f1,f2,u,v,nx,ny,bx,by,1.0f,1.0f,
					  alpha,epsilon_d,diffusivity_type_d,num_iterations_inner,
					  num_iterations_outer,omega,warp_eta,max_warp_levels);
	

	/* read u and v from txt file */
	
	/*
	int i,j;
	FILE *data;
	data = fopen("u.txt","r");
	for ( j = 0;j < ny; j++)
		for (i = 0; i < nx; i++)
			fscanf(data,"%f",&u[i][j]);
	fclose(data);
	data = fopen("v.txt","r");
	for ( j = 0;j < ny; j++)
		for (i = 0; i < nx; i++)
			fscanf(data,"%f",&v[i][j]);
	fclose(data);
	*/
	


	float ref_d,cal_d;
	calculate_errors_2d(u_truth,v_truth,u,v,nx,ny,bx,by,&aae,&al2e,&ref_d,&cal_d);
					  
	compose_image();
	
	// write new computation result to rotating sequence of files
	str = (char*) calloc(2048,sizeof(char));
	writeParamsToStr(str);
	write_ppm_image_comment(get_lastname(),nx,ny,bx,by,of_rgb,1,"%s",str);
	free(str);
	
	// display
	glutPostRedisplay();
}

/*-------------------------------------------------------------------------------------*/

void handleKeyboardspecial(int key, int x, int y)
{   
	change_valf(key,active_param,'a',&alpha);
	change_vali(key,active_param,'d',&diffusivity_type_d);
	change_valf(key,active_param,'e',&warp_eta);
	change_valf(key,active_param,'E',&epsilon_d);
	change_vali(key,active_param,'i',&num_iterations_inner);
	change_vali(key,active_param,'I',&num_iterations_outer);
	change_valf(key,active_param,'l',&max_displacement);
	change_vali(key,active_param,'m',&max_warp_levels);
	change_valf(key,active_param,'O',&omega);
	clamp(0,&num_iterations_inner,10000000);
	clamp(0,&num_iterations_outer,10000000);
	clampf(0,&omega,2);
	clampf(0,&warp_eta,1);
	clamp(0,&max_warp_levels,10000000);
	showParams();    
}

/*-------------------------------------------------------------------------------------*/

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

/*-------------------------------------------------------------------------------------*/

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

/*-------------------------------------------------------------------------------------*/

void ReSizeGLScene(int Width, int Height)
{
	glViewport(0,0,Width,Height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho (0, num_images_side_by_side*nx, ny, 0, -1.0f, 1.0f);
	glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

/*-------------------------------------------------------------------------------------*/

int main (int argc, char* argv[])
{
	int maxgv;
	char * console_string = (char*)calloc(3000,sizeof(char));
	lastclickx=lastclicky=-1;
	active_param='\0';
	num_images_side_by_side=2;
	int no_frontend;
	parse_arg_int   (argc,argv,"bordersizex",&bx,2,console_string);
	parse_arg_int   (argc,argv,"bordersizey",&by,2,console_string);
	parse_arg_float (argc,argv,"alpha",&alpha,100,console_string);
	parse_arg_float (argc,argv,"epsilon_d",&epsilon_d,0.01f,console_string);
	parse_arg_int   (argc,argv,"diffusivity_type_d",&diffusivity_type_d,0,console_string);
	parse_arg_int   (argc,argv,"iterations_inner",&num_iterations_inner,30,console_string);
	parse_arg_int   (argc,argv,"iterations_outer",&num_iterations_outer,5,console_string);
	parse_arg_float (argc,argv,"omega",&omega,1.9f,console_string);
	parse_arg_float (argc,argv,"eta",&warp_eta,0.5f,console_string);
	parse_arg_int   (argc,argv,"max_warp_levels",&max_warp_levels,200,console_string);
	parse_arg_float (argc,argv,"max_displacement",&max_displacement,5.0f,console_string);
	parse_arg_string(argc,argv,"basename",basename,"../images/yos",console_string);
	parse_arg_int   (argc,argv,"no_frontend",&no_frontend,0,console_string);
	char filename[200];
	printf("------\n%s---------\n",console_string);
	
	sprintf(filename,"%s1.pgm",basename);
	f1 = read_pgm_image(filename,&nx,&ny,bx,by,&maxgv);
	
	sprintf(filename,"%s2.pgm",basename);
	f2 = read_pgm_image(filename,&nx,&ny,bx,by,&maxgv);
        
	//printf("f1[100][150] =%f\n",f1[100][150]);
	//printf("f1[150][20] =%f\n",f1[150][20]);
	calloc_multi(1,3,sizeof(float),3,nx,ny,0,bx,by,0,0,0,&of_rgb);
	calloc_multi(4,2,sizeof(float),nx,ny,bx,by,0,0,&u,&v,&u_truth,&v_truth);
	calloc_multi(1,3,sizeof(uchar),num_images_side_by_side*nx,ny,
				 3, 0,0,0, 0,0,0, &p6);
//		printf("f1[150][20] =%f\n",f1[150][20]);
	sprintf(filename,"%st.F",basename);
	read_barron_data(filename,u_truth,v_truth,nx,ny,bx,by);//------------------------------uncomment for computing errors.!!!!!!!!!!!!!!!!!!!!

/*	rewrite the barron file to txt file	*/
//	FILE* UT;
//	FILE* VT;
//	int i,j;
//	UT = fopen("utruth.txt","w");
//	VT = fopen("vtruth.txt","w");
//	for(i = bx; i <nx+bx; i++){
//		for(j = by; j < by+ny; j++)
//		{
//			fprintf(UT,"%f\t",u_truth[i][j]);
//			fprintf(VT,"%f\t",v_truth[i][j]);
//		}
//		fprintf(UT,"\n");
//		fprintf(VT,"\n");
//	}
//	fclose(UT);
//	fclose(VT);

	if(!no_frontend)
	{
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
	}
	else
	{
		compute();
		showParams();
	}

	return(0);
}


#endif

