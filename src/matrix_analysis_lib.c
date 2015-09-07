#ifndef MATRIX_ANALYSIS
#define MATRIX_ANALYSIS

/*---------------------------------------------------------------------------*/
/*             copyright 2011 by Oliver Demetz                               */      
/*                               Mathematical Image Analysis Group           */
/*                               Saarland University                         */
/*---------------------------------------------------------------------------*/


//=================================================

float sum
(
	float **u,
	int nx,
	int ny,
	int bx,
	int by
)
{
	int i,j;
	float s = 0.0f;
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
			s += (u[i][j]);
	return s;
}

//=================================================

float mean_val
(
	float **u,
	int nx,
	int ny,
	int bx,
	int by
)
{
	return sum(u,nx,ny,bx,by) / ((float)(nx*ny));
}

//=================================================

float min_val
(
	float **u,
	int nx,
	int ny,
	int bx,
	int by
)
{
	int i,j;
	float min = u[bx][by];
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
			if(u[i][j]<min)
				min=u[i][j];
	return min;
}

//=================================================

float max_val
(
	float **u,
	int nx,
	int ny,
	int bx,
	int by
)
{
	int i,j;
	float max = u[bx][by];
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
			if(u[i][j]>max)
				max=u[i][j];
	return max;
}

//=================================================

float l2norm
(
	float **u,
	int nx,
	int ny,
	int bx,
	int by
)
{
	int i,j;
	float norm = 0;
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
			norm += u[i][j]*u[i][j];
	return sqrtf(norm);
}

//=================================================

float l1norm
(
	float **u,
	int nx,
	int ny,
	int bx,
	int by
)
{
	int i,j;
	float norm = 0;
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
			norm += fabs(u[i][j]);
	return norm;
}

//=================================================

float var
(
	float **u,
	int nx,
	int ny,
	int bx,
	int by
)
{
	int i,j;
	float mean;
	float var = 0;
	mean = mean_val(u,nx,ny,bx,by);
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
			var += (u[i][j]-mean)*(u[i][j]-mean);
	return var / ((float)(nx*ny));
}


float var_twin
(
	float **u,
	float **v,
	int nx,
	int ny,
	int bx,
	int by
)
{
	int i,j;
	double mean;
	double var = 0;
	mean = (double)(0.5f*(mean_val(u,nx,ny,bx,by)+mean_val(v,nx,ny,bx,by)));
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
			var += (double)( (u[i][j]-mean)*(u[i][j]-mean) 
			               + (v[i][j]-mean)*(v[i][j]-mean) );
	return (float)(var / ((double)(2*nx*ny)));
}

float var_twin_A
(
	float **u,
	float **v,
	float **A11,
	float **A12,
	float **A22,
	int nx,
	int ny,
	int bx,
	int by
)
{
	int i,j;
	float meanu,meanv;
	double var = 0;
	float um,vm;
	meanu = mean_val(u,nx,ny,bx,by);
	meanv = mean_val(v,nx,ny,bx,by);
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
		{
			um = u[i][j]-meanu;
			vm = v[i][j]-meanv;
			var += (double)
				(um*(A11[i][j]*um+A12[i][j]*vm)
				+vm*(A12[i][j]*um+A22[i][j]*vm));
		}
	return (float)(var / ((double)(2*nx*ny)));
}

float cov
(
	float **u,
	float **v,
	int nx,
	int ny,
	int bx,
	int by
)
{
	int i,j;
	float meanu,meanv;
	float cov = 0;
	meanu=mean_val(u,nx,ny,bx,by);
	meanv=mean_val(v,nx,ny,bx,by);
	
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
			cov += (u[i][j]-meanu)*(v[i][j]-meanv);
	return cov/((float)(nx*ny-1));
}

float corr
(
	float **u,
	float **v,
	int nx,
	int ny,
	int bx,
	int by
)
{
	return cov(u,v,nx,ny,bx,by) / 
		sqrtf(var(u,nx,ny,bx,by)*var(v,nx,ny,bx,by));
}

float cov_twin_A
(
	float **u1,
	float **u2,
	float **v1,
	float **v2,
	float **A11,
	float **A12,
	float **A22,
	int nx,
	int ny,
	int bx,
	int by
)
{
	int i,j;
	float meanu1,meanv1,meanu2,meanv2;
	float u1m,u2m,v1m,v2m;
	float cov = 0;
	meanu1=mean_val(u1,nx,ny,bx,by);
	meanu2=mean_val(u2,nx,ny,bx,by);
	meanv1=mean_val(v1,nx,ny,bx,by);
	meanv2=mean_val(v2,nx,ny,bx,by);
	
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
		{
			u1m=u1[i][j]-meanu1; 
			u2m=u2[i][j]-meanu2;
			v1m=v1[i][j]-meanv1;
			v2m=v2[i][j]-meanv2;
			cov += u1m*(v1m*A11[i][j]+v2m*A12[i][j]) 
			          + u2m*(v1m*A12[i][j]+v2m*A22[i][j]);
		}
	return cov/((float)(2*nx*ny-1));
}

float corr_twin_A
(
	float **u1,
	float **u2,
	float **v1,
	float **v2,
	float **A11,
	float **A12,
	float **A22,
	int nx,
	int ny,
	int bx,
	int by
)
{
	return cov_twin_A(u1,u2,v1,v2,A11,A12,A22,nx,ny,bx,by) /
		sqrtf( var_twin_A(u1,u2,A11,A12,A22,nx,ny,by,by)
		      *var_twin_A(v1,v2,A11,A12,A22,nx,ny,bx,by) );
}


#endif

