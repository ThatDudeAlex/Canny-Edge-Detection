#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define  PICSIZE 256
#define  MAXMASK 100

int    pic[PICSIZE][PICSIZE];
double outpicx[PICSIZE][PICSIZE];
double outpicy[PICSIZE][PICSIZE];
double xmask[MAXMASK][MAXMASK];
double ymask[MAXMASK][MAXMASK];
int histogram[PICSIZE];
double magPic[PICSIZE][PICSIZE];
double peaksPic[PICSIZE][PICSIZE];
double finalPic[PICSIZE][PICSIZE];


int main(int argc, char **argv)
{
  int     i,j,p,q, mr,centx, centy, moretodo;
  double  maskval,maxival,minival,maxval, sig, sum1, sum2, slope;
  double high, low, percent, cutoff = 0, area = 0;
  FILE    *fo1, *fo2, *fo3, *fp1, *fopen();

  fp1=fopen(argv[1],"rb");
  fseek(fp1, 15, SEEK_SET );

  fo1=fopen("cannyMag.pgm","wb");
  fprintf(fo1,"P5\n256 256\n255\n");

  fo2=fopen("cannyPeaks.pgm","wb");
  fprintf(fo2,"P5\n256 256\n255\n");

  fo3=fopen("cannyFinal.pgm","wb");
  fprintf(fo3,"P5\n256 256\n255\n");

   sig = atof(argv[2]);

   percent = atof(argv[3]) / 100;

   mr = (int)(sig * 3);
   centx = (MAXMASK / 2);
   centy = (MAXMASK / 2);

   fseek(fp1, 15, SEEK_SET );


  // read pic
  for (i=0;i<256;i++)
  { for (j=0;j<256;j++)
          {
            pic[i][j]  =  getc (fp1);
          }
  }

  //-------------part 1

  // make Xmask and Ymask
   for (p=-mr;p<=mr;p++)
   {  for (q=-mr;q<=mr;q++)
      {
         maskval = (q)*(exp(-1*(((p*p)+(q*q))/(2*(sig*sig)))));
         xmask[p+centy][q+centx] = maskval;

         maskval = (p)*(exp(-1*(((p*p)+(q*q))/(2*(sig*sig)))));
         ymask[p+centy][q+centx] = maskval;
      }
   }

   // convolution to get Xpic and Ypic
   for (i=mr;i<=255-mr;i++)
   {
     for (j=mr;j<=255-mr;j++)
     {
        sum1 = 0;
        sum2 = 0;
        for (p=-mr;p<=mr;p++)
        {
           for (q=-mr;q<=mr;q++)
           {
              sum1 += pic[i+p][j+q] * xmask[p+centy][q+centx];
              sum2 += pic[i+p][j+q] * ymask[p+centy][q+centx];
           }
        }
        outpicx[i][j] = sum1;
        outpicy[i][j] = sum2;
     }
   }

   // compute magnitide
   maxival = 0;
   for (i=mr;i<256-mr;i++)
   { for (j=mr;j<256-mr;j++)
     {
        magPic[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                 (outpicy[i][j]*outpicy[i][j])));
        if (magPic[i][j] > maxival)
           maxival = magPic[i][j];

      }
   }

   // print to file magnitude image
    for (i=0;i<256;i++)
    {
       for (j=0;j<256;j++)
       {
          magPic[i][j] = (magPic[i][j] / maxival) * 255;
          fprintf(fo1,"%c",(char)((int)(magPic[i][j])));
       }
    }

   // part 2

   // checks pairs of 3 looking for peaks
   for(i=mr;i<256-mr;i++)
   {
     for(j=mr;j<256-mr;j++)
     {

       if((outpicx[i][j]) == 0.0)
       {
         outpicx[i][j] = .00001;
       }

      slope = outpicy[i][j]/outpicx[i][j];

      if( (slope <= .4142)&&(slope > -.4142))
      {
         if((magPic[i][j] > magPic[i][j-1])&&(magPic[i][j] > magPic[i][j+1]))
         {
            peaksPic[i][j] = 255;
         }
      }
      else if( (slope <= 2.4142)&&(slope > .4142))
      {
         if((magPic[i][j] > magPic[i-1][j-1])&&(magPic[i][j] > magPic[i+1][j+1]))
         {
             peaksPic[i][j] = 255;
         }
      }
      else if( (slope <= -.4142)&&(slope > -2.4142))
      {
         if((magPic[i][j] > magPic[i+1][j-1])&&(magPic[i][j] > magPic[i-1][j+1]))
         {
             peaksPic[i][j] = 255;
         }
       }
       else
       {
         if((magPic[i][j] > magPic[i-1][j])&&(magPic[i][j] > magPic[i+1][j]))
         {
             peaksPic[i][j] = 255;
         }
       }
     }
   }

   // print to file peaks image
    for (i=0;i<256;i++)
    {
       for (j=0;j<256;j++)
       {
          fprintf(fo2,"%c",(char)((int)(peaksPic[i][j])));
       }
    }

    // part 4--------------------------------

    for(i = 0; i < 256; i++)
    {
      for(j = 0; j< 256; j ++)
      {
        histogram[(int)magPic[i][j]]++;
      }
    }

    cutoff = PICSIZE * PICSIZE * percent;

    for(i = 255; i >= 1; i--)
    {
      area += histogram[i];

      if(area > cutoff)
      {
        break;
      }

    }

    high = i;
    low = high * 0.35;

    printf("High Threshold: %.2f\nLow Threshold: %.2f", high, low);


    // part 3---------------

    for(i = 0; i < 256; i++)
    {
      for(j = 0; j < 256; j++)
      {
        if(peaksPic[i][j] == 255)
        {
          if(magPic[i][j] > high)
          {
            peaksPic[i][j] = 0;
            finalPic[i][j] = 255;
          }
          else if(magPic[i][j] < low)
          {
            peaksPic[i][j] = finalPic[i][j] = 0;
          }
        }
      }
    }


moretodo = 1;
while(moretodo == 1)
{
  moretodo = 0;

  for(i = 0; i < 256; i++)
  {
    for(j = 0; j < 256; j++)
    {
      if(peaksPic[i][j] == 255)
      {
        for(p = -1; p <= 1; p++)
        {
          for(q = -1; q <= 1; q++)
          {
            if(finalPic[i + p][j + q] == 255)
            {
              peaksPic[i][j] = 0;
              finalPic[i][j] = 255;
              moretodo = 1;
            }
          }
        }
      }
    }
  }
}


    for (i=0;i<256;i++)
    {
       for (j=0;j<256;j++)
       {
          fprintf(fo3,"%c",(char)((int)(finalPic[i][j])));
       }
    }


    fclose(fo1);
    fclose(fo2);
    fclose(fo3);
    fclose(fp1);


}
