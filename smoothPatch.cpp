#include "math.h"
#include <cstdio>
#include <cstring>
#include <cstdint>


using namespace std;

void smoothPatchC(float* vertices, int VertexN,int64_t* faces, int FacesN);

void smoothPatchC(float* vertices, int VertexN,int64_t* faces, int FacesN){
    /* All inputs */
    double *FacesA, *FacesB, *FacesC;
    double* VerticesX;
    double *VerticesY;
    double *VerticesZ;
    int Iterations, Lambda;
    /* All outputs */
    double *VerticesOX, *VerticesOY, *VerticesOZ;

    /* All outputs, Vertex storage*/
    double *VerticesNX, *VerticesNY, *VerticesNZ;
    double *VerticesN2X, *VerticesN2Y, *VerticesN2Z;
    
        
    /* Temporary Weights */
    double *VerticesW ;
      
    
    /* Point Update temporary storage */
    double Ux, Uy, Uz;
    
    /* 1D Index  */
    int index0, index1, index2;
            
    /* Edge coordinates and lenght */
    double e0x, e0y, e0z, e0l;
    double e1x, e1y, e1z, e1l;
    double e2x, e2y, e2z, e2l;
    
    /* Swap point variable */
    double *t;
    int swap=0;
    
    /* Loop variable */
    int i,j;

   
   VerticesX = (double *)malloc( VertexN/3 * sizeof(double) );
   VerticesY = (double *)malloc( VertexN/3 * sizeof(double) );
   VerticesZ = (double *)malloc( VertexN/3 * sizeof(double) );
   FacesA = (double *)malloc( FacesN/4 * sizeof(double) );
   FacesB = (double *)malloc( FacesN/4 * sizeof(double) );
   FacesC = (double *)malloc( FacesN/4 * sizeof(double) );

   /* Read all inputs (faces and vertices) */

    for(int i = 0; i < VertexN; i+=3){
       
        VerticesX[(int)(i/3)] = vertices[i];
        VerticesY[(int)(i/3)] = vertices[i+1];
        VerticesZ[(int)(i/3)] = vertices[i+2];
    }

    for(int i = 0; i < FacesN; i+=4){

        FacesA[(int)(i/4)] = faces[i+1];
        FacesB[(int)(i/4)] = faces[i+2];
        FacesC[(int)(i/4)] = faces[i+3];
        
    }

   Iterations=20;
   Lambda=1;

   /* Intern vertices storage */
   VerticesW = (double *)malloc( VertexN/3 * sizeof(double) );
   VerticesNX = (double *)malloc( VertexN/3 * sizeof(double) );
   VerticesNY = (double *)malloc( VertexN/3 * sizeof(double) );
   VerticesNZ = (double *)malloc( VertexN/3 * sizeof(double) );
   VerticesN2X = (double *)malloc( VertexN/3 * sizeof(double) );
   VerticesN2Y = (double *)malloc( VertexN/3 * sizeof(double) );
   VerticesN2Z = (double *)malloc( VertexN/3 * sizeof(double) );
   
   /* Copy input arrays to ouput vertice arrays */
   memcpy( VerticesNX,VerticesX,VertexN/3 * sizeof(double));
   memcpy( VerticesNY,VerticesY,VertexN/3 * sizeof(double));
   memcpy( VerticesNZ,VerticesZ,VertexN/3 * sizeof(double));
   

   for (j=0; j<Iterations; j++)
   {
       /* Clean the weights */
       for (i=0; i<VertexN/3; i++) { VerticesW[i] =0;  VerticesN2X[i]=0; VerticesN2Y[i]=0; VerticesN2Z[i]=0; }

       /* Calculate all face normals and angles */

       #pragma parallel for
       for (i=0; i<FacesN/4; i++)
       {
           /* Get indices of face vertices */
           index0=(int)FacesA[i];
           index1=(int)FacesB[i];
           index2=(int)FacesC[i];

           /* Calculate edge lengths */
           e0x=VerticesNX[index0]-VerticesNX[index1];  
           e0y=VerticesNY[index0]-VerticesNY[index1];  
           e0z=VerticesNZ[index0]-VerticesNZ[index1];
           e1x=VerticesNX[index1]-VerticesNX[index2];  
           e1y=VerticesNY[index1]-VerticesNY[index2];  
           e1z=VerticesNZ[index1]-VerticesNZ[index2];
           e2x=VerticesNX[index2]-VerticesNX[index0];  
           e2y=VerticesNY[index2]-VerticesNY[index0];  
           e2z=VerticesNZ[index2]-VerticesNZ[index0];
           e0l=1 / (sqrt(e0x*e0x + e0y*e0y + e0z*e0z)+Lambda);
           e1l=1 / (sqrt(e1x*e1x + e1y*e1y + e1z*e1z)+Lambda);
           e2l=1 / (sqrt(e2x*e2x + e2y*e2y + e2z*e2z)+Lambda);

           VerticesN2X[index0]+=VerticesNX[index1]*e0l;
           VerticesN2Y[index0]+=VerticesNY[index1]*e0l;
           VerticesN2Z[index0]+=VerticesNZ[index1]*e0l;
           VerticesW[index0]+=e0l;

           VerticesN2X[index1]+=VerticesNX[index0]*e0l; 
           VerticesN2Y[index1]+=VerticesNY[index0]*e0l;
           VerticesN2Z[index1]+=VerticesNZ[index0]*e0l;
           VerticesW[index1]+=e0l;


           VerticesN2X[index1]+=VerticesNX[index2]*e1l;
           VerticesN2Y[index1]+=VerticesNY[index2]*e1l;
           VerticesN2Z[index1]+=VerticesNZ[index2]*e1l;
           VerticesW[index1]+=e1l;

           VerticesN2X[index2]+=VerticesNX[index1]*e1l;
           VerticesN2Y[index2]+=VerticesNY[index1]*e1l;
           VerticesN2Z[index2]+=VerticesNZ[index1]*e1l;
           VerticesW[index2]+=e1l;

           VerticesN2X[index2]+=VerticesNX[index0]*e2l;
           VerticesN2Y[index2]+=VerticesNY[index0]*e2l; 
           VerticesN2Z[index2]+=VerticesNZ[index0]*e2l;
           VerticesW[index2]+=e2l;

           VerticesN2X[index0]+=VerticesNX[index2]*e2l; 
           VerticesN2Y[index0]+=VerticesNY[index2]*e2l; 
           VerticesN2Z[index0]+=VerticesNZ[index2]*e2l;
           VerticesW[index0]+=e2l;
       }

       /* Normalize the Vertices */
       #pragma parallel for
       for (i=0; i<VertexN/3; i++)
       {
           Ux=0; Uy=0; Uz=0;
           Ux=VerticesN2X[i]/VerticesW[i];
           Uy=VerticesN2Y[i]/VerticesW[i];
           Uz=VerticesN2Z[i]/VerticesW[i];

           Ux=Ux-VerticesNX[i]; 
           Uy=Uy-VerticesNY[i]; 
           Uz=Uz-VerticesNZ[i]; 

           VerticesN2X[i]=VerticesNX[i]+Ux*Lambda; 
           VerticesN2Y[i]=VerticesNY[i]+Uy*Lambda; 
           VerticesN2Z[i]=VerticesNZ[i]+Uz*Lambda;
       }
       
       /* Swap the variables */
       t=VerticesNX; VerticesNX=VerticesN2X; VerticesN2X=t;
       t=VerticesNY; VerticesNY=VerticesN2Y; VerticesN2Y=t;
       t=VerticesNZ; VerticesNZ=VerticesN2Z; VerticesN2Z=t;

       /* Swap output variable */
       if(swap==0) { swap=1; } else { swap=0; }
   }

   if(swap==0)
   {


        for(int i = 0; i < VertexN; i+=3){

            vertices[i] = VerticesN2X[(int)(i/3)];
            vertices[i+1] = VerticesN2Y[(int)(i/3)];
            vertices[i+2] = VerticesN2Z[(int)(i/3)];
        }

   }
   else
   {
       
        for(int i = 0; i < VertexN; i+=3){

            vertices[i] = VerticesNX[(int)(i/3)];
            vertices[i+1] = VerticesNY[(int)(i/3)];
            vertices[i+2] = VerticesNZ[(int)(i/3)];
        }
   }
           
   
}

extern "C" {
    void smoothPatch(float* vertices, int VertexN, int64_t* faces, int FacesN)
    {
        return smoothPatchC(vertices,VertexN,faces,FacesN);
    }
}
 


