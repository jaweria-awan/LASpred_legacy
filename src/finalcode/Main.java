/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package finalcode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author Javeria
 */
public class Main {
    static final int ROWS = 20;
    static final int COLS = 21;
    static final int  [][] a2 = new  int  [ROWS][COLS];
    static final double    [][] af = new  double  [ROWS][COLS];
 static    int count=0;


    /**
     * @param args the command line arguments
     */
      public static void main(String[] args) throws FileNotFoundException, IOException {
        BufferedReader reader = new BufferedReader(new FileReader("E:/e-coli15pos.txt"));//give read file path//REPLACE WITH INPUT FROM GUI
       BufferedWriter writer = new BufferedWriter(new FileWriter("E:/outpute-coli.txt"));//give write file path
       String stri="";
       String line="";
       double Threshold=212.368;
        PSSMmaker();
        freqcal(a2);
      line=reader.readLine();
    while(line!=null)
    {
        
        if(line.contains(">")){
            stri="";
        writer.write(line);
        writer.newLine();
        line=reader.readLine();
        if(line==null)
            break;
        while (!(line.contains(">"))){
            stri=stri+""+line;
            line=reader.readLine();
            if(line==null)
                break;
            
        }
        }
           writer.write(stri);
           writer.newLine();
       
          
            for(int i=0;i<stri.length();i++)
                
              {
               
               if (stri.charAt(i)=='K')
               {
               int pos=i;
               String posi=Integer.toString(pos);
               int forwerd=pos+11;
               int backward=pos-10;
               if(backward>0&&forwerd<stri.length())
               {
                   writer.write("Position Of K");
                   writer.write(posi+"      ");
                   //writer.newLine();
                       String substri=stri.substring(backward,forwerd) ;
                        writer.write(substri+"      ");
                       // writer.newLine();
                        double AbsEntp=0.0;
                        double AvgVolRes=0.0;
                        double AvgNonBondEng=0.0;
                        double Heatcap=0.0;
                        double IsoElectricPoint=0.0;
                        double Seqfreq=0.0;
                        double Chufasman=0.0;
                        double PSSM=0.0;
                        double Dimer=0.0;
                        double IntgScore=0.0;
                       AbsEntp=AbsolutEentropy (substri);
                       AvgVolRes= AverageVolumesOfResidues (substri);
                       AvgNonBondEng= AverageNonBondedEnergyPerResidue (substri);
                       Heatcap= HeatCapacity (substri);
                       IsoElectricPoint=  IsoelectricPoint  (substri);
                       Seqfreq= SequenceFrequency (substri);
                       Chufasman= TheChouFasmanParameterOfTheCoilConformation  (substri);
                       PSSM=scoring(substri);
                       Dimer= DimerScore (substri);
                        IntgScore=0.0;
                       IntgScore=((AbsEntp)+(AvgVolRes)+(AvgNonBondEng)+(Heatcap)+(IsoElectricPoint)+(Seqfreq)+(Chufasman)+(PSSM)+(Dimer));
                       
                       if(IntgScore>Threshold){
                           String intgscore=Double.toString(IntgScore);
                           writer.write(intgscore);
                           writer.newLine();
                           //writer.write("The site is a potential acetylated Site");
                           //writer.newLine();
                           
                       }
                      /* else{
                           writer.write("The site is not a potential acetylated Site");
                           writer.newLine();
                       }*/
                      
                }
            }
            }
          }

    
       writer.close();
        reader.close();

}

  
        // TODO code application logic here
    
     public static double AbsolutEentropy (String stri) throws FileNotFoundException, IOException
    {
      
       
       double Ascore=30.88;
       double Lscore=50.62;
       double Rscore= 68.43;
       double Kscore=63.21;
       double Nscore=41.70;
       double Mscore=55.32;
       double Dscore=40.66;
       double Fscore=51.06;
       double Cscore=53.83;
       double Pscore=39.21;
       double Qscore=46.62;
       double Sscore=35.65;
       double Escore=44.98;
       double Tscore=36.50;
       double Gscore=24.74;
       double Wscore=60.00;
       double Hscore=65.99;
       double Yscore=51.15;
       double Iscore=49.71;
       double Vscore=42.75;
       double  fscore=0.0;
       double coefabsentp =-0.026;
    
           for(int i=0;i<stri.length();i++){
              
                            if( stri.charAt(i)=='A')
                              fscore+=Ascore;
                            else   if( stri.charAt(i)=='G')
                              fscore+=Gscore;
                            else  if( stri.charAt(i)=='V')
                              fscore+=Vscore;
                            else if( stri.charAt(i)=='L')
                              fscore+=Lscore;
                            else if( stri.charAt(i)=='I')
                             fscore+=Iscore;
                             else if( stri.charAt(i)=='S')
                             fscore+=Sscore;
                             else if( stri.charAt(i)=='T')
                             fscore+=Tscore;
                             else if( stri.charAt(i)=='C')
                            fscore+=Cscore;
                             else if( stri.charAt(i)=='M')
                            fscore+=Mscore;
                             else  if( stri.charAt(i)=='P')
                             fscore+=Pscore;
                             else   if( stri.charAt(i)=='H')
                             fscore+=Hscore;
                             else if( stri.charAt(i)=='R')
                              fscore+=Rscore;
                             else  if( stri.charAt(i)=='N')
                            fscore+=Nscore;
                             else if( stri.charAt(i)=='Q')
                             fscore+=Qscore;
                             else  if( stri.charAt(i)=='E')
                             fscore+=Escore;
                             else if( stri.charAt(i)=='D')
                            fscore+=Dscore;
                             else  if( stri.charAt(i)=='F')
                            fscore+=Fscore;
                             else  if( stri.charAt(i)=='W')
                            fscore+=Wscore;
                             else  if( stri.charAt(i)=='Y')
                             fscore+=Yscore;
                             else if ( stri.charAt(i)=='K')
                             fscore+=Kscore;
          }
           
            fscore=(fscore*coefabsentp );
             fscore=roundtwofour(fscore);
            return fscore;
           
             //String line=Double.toString(fscore);
             
          
       
   
     }
     public static double  AverageVolumesOfResidues (String stri) throws 
FileNotFoundException, IOException
    {
       double Ascore=91.5;
       double Lscore=163.4;
       double Rscore=196.1;
       double Kscore=162.5;
       double Nscore=138.3;
       double Mscore=165.9;
       double Dscore=135.2;
       double Fscore=198.8;
       double Cscore=114.4;
       double Pscore=123.4;
       double Qscore=156.4;
       double Sscore=102.0;
       double Escore=154.6;
       double Tscore=126.0;
       double Gscore=67.5;
       double Wscore=209.8;
       double Hscore=163.2;
       double Yscore=237.2;
       double Iscore=162.6;
       double Vscore=138.4;
       double  fscore=0.0;
       double coefavgvol=(-0.009);
    
           for(int i=0;i<stri.length();i++){
              
                            if( stri.charAt(i)=='A')
                              fscore+=Ascore;
                            else   if( stri.charAt(i)=='G')
                              fscore+=Gscore;
                            else  if( stri.charAt(i)=='V')
                              fscore+=Vscore;
                            else if( stri.charAt(i)=='L')
                              fscore+=Lscore;
                            else if( stri.charAt(i)=='I')
                             fscore+=Iscore;
                             else if( stri.charAt(i)=='S')
                             fscore+=Sscore;
                             else if( stri.charAt(i)=='T')
                             fscore+=Tscore;
                             else if( stri.charAt(i)=='C')
                            fscore+=Cscore;
                             else if( stri.charAt(i)=='M')
                            fscore+=Mscore;
                             else  if( stri.charAt(i)=='P')
                             fscore+=Pscore;
                             else   if( stri.charAt(i)=='H')
                             fscore+=Hscore;
                             else if( stri.charAt(i)=='R')
                              fscore+=Rscore;
                             else  if( stri.charAt(i)=='N')
                            fscore+=Nscore;
                             else if( stri.charAt(i)=='Q')
                             fscore+=Qscore;
                             else  if( stri.charAt(i)=='E')
                             fscore+=Escore;
                             else if( stri.charAt(i)=='D')
                            fscore+=Dscore;
                             else  if( stri.charAt(i)=='F')
                            fscore+=Fscore;
                             else  if( stri.charAt(i)=='W')
                            fscore+=Wscore;
                             else  if( stri.charAt(i)=='Y')
                             fscore+=Yscore;
                             else if ( stri.charAt(i)=='K')
                             fscore+=Kscore;
          }
          
            fscore=(fscore*coefavgvol);
          
            fscore=roundtwofour(fscore);
    
            return fscore;
     }
       public static double AverageNonBondedEnergyPerResidue (String stri) throws 
FileNotFoundException, IOException
    {
       double Ascore=-7.020 ;
       double Lscore=-10.520;
       double Rscore=-10.131;
       double Kscore=-9.666;
       double Nscore=-9.424;
       double Mscore=-10.424;
       double Dscore=-9.296;
       double Fscore=-12.485;
       double Cscore=-8.190;
       double Pscore=-8.652;
       double Qscore=-10.044;
       double Sscore=-7.782;
       double Escore=-10.467;
       double Tscore=-8.764;
       double Gscore=-5.456;
       double Wscore=-14.420;
       double Hscore=-12.150;
       double Yscore=-12.360;
       double Iscore=-9.512;
       double Vscore=-8.778;
       double  fscore=0.0;
       double coefavgnon= -0.695;
    
           for(int i=0;i<stri.length();i++){
              
                            if( stri.charAt(i)=='A')
                              fscore+=Ascore;
                            else   if( stri.charAt(i)=='G')
                              fscore+=Gscore;
                            else  if( stri.charAt(i)=='V')
                              fscore+=Vscore;
                            else if( stri.charAt(i)=='L')
                              fscore+=Lscore;
                            else if( stri.charAt(i)=='I')
                             fscore+=Iscore;
                             else if( stri.charAt(i)=='S')
                             fscore+=Sscore;
                             else if( stri.charAt(i)=='T')
                             fscore+=Tscore;
                             else if( stri.charAt(i)=='C')
                            fscore+=Cscore;
                             else if( stri.charAt(i)=='M')
                            fscore+=Mscore;
                             else  if( stri.charAt(i)=='P')
                             fscore+=Pscore;
                             else   if( stri.charAt(i)=='H')
                             fscore+=Hscore;
                             else if( stri.charAt(i)=='R')
                              fscore+=Rscore;
                             else  if( stri.charAt(i)=='N')
                            fscore+=Nscore;
                             else if( stri.charAt(i)=='Q')
                             fscore+=Qscore;
                             else  if( stri.charAt(i)=='E')
                             fscore+=Escore;
                             else if( stri.charAt(i)=='D')
                            fscore+=Dscore;
                             else  if( stri.charAt(i)=='F')
                            fscore+=Fscore;
                             else  if( stri.charAt(i)=='W')
                            fscore+=Wscore;
                             else  if( stri.charAt(i)=='Y')
                             fscore+=Yscore;
                             else if ( stri.charAt(i)=='K')
                             fscore+=Kscore;
          }
           
            fscore=(fscore*coefavgnon);
             fscore=roundtwofour(fscore);
            return fscore;
       }
     
       public static double  HeatCapacity (String stri) throws 
FileNotFoundException, IOException
    {
       
       double Ascore=29.22;
       double Lscore=48.03;
       double Rscore=26.37;
       double Kscore=57.10;
       double Nscore=38.30;
       double Mscore=69.32;
       double Dscore=37.09;
       double Fscore=48.52;
       double Cscore=50.70;
       double Pscore=36.13;
       double Qscore=44.02;
       double Sscore=32.40;
       double Escore=41.84;
       double Tscore=35.20;
       double Gscore=23.71;
       double Wscore=56.92;
       double Hscore=59.64;
       double Yscore=51.73;
       double Iscore=45.00;
       double Vscore=40.35;
       double  fscore=0.0;
       double coefHC= 0.019;
    
           for(int i=0;i<stri.length();i++){
              
                            if( stri.charAt(i)=='A')
                              fscore+=Ascore;
                            else   if( stri.charAt(i)=='G')
                              fscore+=Gscore;
                            else  if( stri.charAt(i)=='V')
                              fscore+=Vscore;
                            else if( stri.charAt(i)=='L')
                              fscore+=Lscore;
                            else if( stri.charAt(i)=='I')
                             fscore+=Iscore;
                             else if( stri.charAt(i)=='S')
                             fscore+=Sscore;
                             else if( stri.charAt(i)=='T')
                             fscore+=Tscore;
                             else if( stri.charAt(i)=='C')
                            fscore+=Cscore;
                             else if( stri.charAt(i)=='M')
                            fscore+=Mscore;
                             else  if( stri.charAt(i)=='P')
                             fscore+=Pscore;
                             else   if( stri.charAt(i)=='H')
                             fscore+=Hscore;
                             else if( stri.charAt(i)=='R')
                              fscore+=Rscore;
                             else  if( stri.charAt(i)=='N')
                            fscore+=Nscore;
                             else if( stri.charAt(i)=='Q')
                             fscore+=Qscore;
                             else  if( stri.charAt(i)=='E')
                             fscore+=Escore;
                             else if( stri.charAt(i)=='D')
                            fscore+=Dscore;
                             else  if( stri.charAt(i)=='F')
                            fscore+=Fscore;
                             else  if( stri.charAt(i)=='W')
                            fscore+=Wscore;
                             else  if( stri.charAt(i)=='Y')
                             fscore+=Yscore;
                             else if ( stri.charAt(i)=='K')
                             fscore+=Kscore;
          }
            
            fscore=(fscore*coefHC);
            fscore=roundtwofour(fscore);
            return fscore;
       }
       public static double  IsoelectricPoint  (String stri) throws 
FileNotFoundException, IOException
    {
       
       double Ascore=6.00;
       double Lscore=5.98;
       double Rscore=10.76;
       double Kscore=9.74;
       double Nscore=5.41;
       double Mscore=5.74;
       double Dscore=2.77;
       double Fscore=5.48;
       double Cscore=5.05;
       double Pscore=6.30;
       double Qscore=5.65;
       double Sscore=5.68;
       double Escore=3.22;
       double Tscore=5.66;
       double Gscore=5.97;
       double Wscore=5.89;
       double Hscore=7.59;
       double Yscore=5.66;
       double Iscore=6.02;
       double Vscore=5.96;
       double  fscore=0.0;
       double iep=0.131;
    
           for(int i=0;i<stri.length();i++){
              
                            if( stri.charAt(i)=='A')
                              fscore+=Ascore;
                            else   if( stri.charAt(i)=='G')
                              fscore+=Gscore;
                            else  if( stri.charAt(i)=='V')
                              fscore+=Vscore;
                            else if( stri.charAt(i)=='L')
                              fscore+=Lscore;
                            else if( stri.charAt(i)=='I')
                             fscore+=Iscore;
                             else if( stri.charAt(i)=='S')
                             fscore+=Sscore;
                             else if( stri.charAt(i)=='T')
                             fscore+=Tscore;
                             else if( stri.charAt(i)=='C')
                            fscore+=Cscore;
                             else if( stri.charAt(i)=='M')
                            fscore+=Mscore;
                             else  if( stri.charAt(i)=='P')
                             fscore+=Pscore;
                             else   if( stri.charAt(i)=='H')
                             fscore+=Hscore;
                             else if( stri.charAt(i)=='R')
                              fscore+=Rscore;
                             else  if( stri.charAt(i)=='N')
                            fscore+=Nscore;
                             else if( stri.charAt(i)=='Q')
                             fscore+=Qscore;
                             else  if( stri.charAt(i)=='E')
                             fscore+=Escore;
                             else if( stri.charAt(i)=='D')
                            fscore+=Dscore;
                             else  if( stri.charAt(i)=='F')
                            fscore+=Fscore;
                             else  if( stri.charAt(i)=='W')
                            fscore+=Wscore;
                             else  if( stri.charAt(i)=='Y')
                             fscore+=Yscore;
                             else if ( stri.charAt(i)=='K')
                             fscore+=Kscore;
          }
            
            fscore=(fscore*iep);
            fscore=roundtwofour(fscore);
            return fscore;
       }
       public static double  SequenceFrequency (String stri) throws 
FileNotFoundException, IOException
    {
       double Ascore=685;
       double Lscore=581;
       double Rscore=382;
       double Kscore=575;
       double Nscore=397;
       double Mscore=132;
       double Dscore=400;
       double Fscore=303;
       double Cscore=241;
       double Pscore=366;
       double Qscore=313;
       double Sscore=593;
       double Escore=427;
       double Tscore=490;
       double Gscore=707;
       double Wscore=99;
       double Hscore=155;
       double Yscore=292;
       double Iscore=394;
       double Vscore=553;
       double  fscore=0.0;
       double sfeq= 0.003;
    
           for(int i=0;i<stri.length();i++){
              
                            if( stri.charAt(i)=='A')
                              fscore+=Ascore;
                            else   if( stri.charAt(i)=='G')
                              fscore+=Gscore;
                            else  if( stri.charAt(i)=='V')
                              fscore+=Vscore;
                            else if( stri.charAt(i)=='L')
                              fscore+=Lscore;
                            else if( stri.charAt(i)=='I')
                             fscore+=Iscore;
                             else if( stri.charAt(i)=='S')
                             fscore+=Sscore;
                             else if( stri.charAt(i)=='T')
                             fscore+=Tscore;
                             else if( stri.charAt(i)=='C')
                            fscore+=Cscore;
                             else if( stri.charAt(i)=='M')
                            fscore+=Mscore;
                             else  if( stri.charAt(i)=='P')
                             fscore+=Pscore;
                             else   if( stri.charAt(i)=='H')
                             fscore+=Hscore;
                             else if( stri.charAt(i)=='R')
                              fscore+=Rscore;
                             else  if( stri.charAt(i)=='N')
                            fscore+=Nscore;
                             else if( stri.charAt(i)=='Q')
                             fscore+=Qscore;
                             else  if( stri.charAt(i)=='E')
                             fscore+=Escore;
                             else if( stri.charAt(i)=='D')
                            fscore+=Dscore;
                             else  if( stri.charAt(i)=='F')
                            fscore+=Fscore;
                             else  if( stri.charAt(i)=='W')
                            fscore+=Wscore;
                             else  if( stri.charAt(i)=='Y')
                             fscore+=Yscore;
                             else if ( stri.charAt(i)=='K')
                             fscore+=Kscore;
          }
           
            fscore=(fscore*sfeq);
             fscore=roundtwofour(fscore);
            return fscore;
       }
       public static double  TheChouFasmanParameterOfTheCoilConformation  (String stri) throws 
FileNotFoundException, IOException
    {
       double Ascore=0.71;
       double Lscore=0.69;
       double Rscore=1.06;
       double Kscore=0.99;
       double Nscore=1.37;
       double Mscore=0.59;
       double Dscore=1.21;
       double Fscore=0.71;
       double Cscore=1.19;
       double Pscore=1.61;
       double Qscore=0.87;
       double Sscore=1.34;
       double Escore=0.84;
       double Tscore=1.08;
       double Gscore=1.52;
       double Wscore=0.76;
       double Hscore=1.07;
       double Yscore=1.07;
       double Iscore=0.66;
       double Vscore=0.63;
       double  fscore=0.0;
       double chu=1.781;
    
           for(int i=0;i<stri.length();i++){
              
                            if( stri.charAt(i)=='A')
                              fscore+=Ascore;
                            else   if( stri.charAt(i)=='G')
                              fscore+=Gscore;
                            else  if( stri.charAt(i)=='V')
                              fscore+=Vscore;
                            else if( stri.charAt(i)=='L')
                              fscore+=Lscore;
                            else if( stri.charAt(i)=='I')
                             fscore+=Iscore;
                             else if( stri.charAt(i)=='S')
                             fscore+=Sscore;
                             else if( stri.charAt(i)=='T')
                             fscore+=Tscore;
                             else if( stri.charAt(i)=='C')
                            fscore+=Cscore;
                             else if( stri.charAt(i)=='M')
                            fscore+=Mscore;
                             else  if( stri.charAt(i)=='P')
                             fscore+=Pscore;
                             else   if( stri.charAt(i)=='H')
                             fscore+=Hscore;
                             else if( stri.charAt(i)=='R')
                              fscore+=Rscore;
                             else  if( stri.charAt(i)=='N')
                            fscore+=Nscore;
                             else if( stri.charAt(i)=='Q')
                             fscore+=Qscore;
                             else  if( stri.charAt(i)=='E')
                             fscore+=Escore;
                             else if( stri.charAt(i)=='D')
                            fscore+=Dscore;
                             else  if( stri.charAt(i)=='F')
                            fscore+=Fscore;
                             else  if( stri.charAt(i)=='W')
                            fscore+=Wscore;
                             else  if( stri.charAt(i)=='Y')
                             fscore+=Yscore;
                             else if ( stri.charAt(i)=='K')
                             fscore+=Kscore;
          }
            
            fscore=(fscore*chu);
            //fscore=roundtwofour(fscore);
            return fscore;
       }
       
       static void PSSMmaker() throws FileNotFoundException, IOException {
      BufferedReader readeri = new BufferedReader(new FileReader("E:/InputFileForHumanPosi.txt"));//give read file path
       
      String s=">";
      String str=null;
        
        while ((str = readeri.readLine()) != null) {
       if  (!str.startsWith(s))
       {
           count ++;
           for(int i=0;i<COLS;i++)
                   {
                             if( str.charAt(i)=='A')
                              a2[0][i]+=1;
                            else   if( str.charAt(i)=='G')
                              a2[1][i]+=1;
                            else  if( str.charAt(i)=='V')
                              a2[2][i]+=1;
                            else if( str.charAt(i)=='L')
                              a2[3][i]+=1;
                            else if( str.charAt(i)=='I')
                             a2[4][i]+=1;
                             else if( str.charAt(i)=='S')
                             a2[5][i]+=1;
                             else if( str.charAt(i)=='T')
                             a2[6][i]+=1;
                             else if( str.charAt(i)=='C')
                            a2[7][i]+=1;
                             else if( str.charAt(i)=='M')
                            a2[8][i]+=1;
                             else  if( str.charAt(i)=='P')
                             a2[9][i]+=1;
                             else   if( str.charAt(i)=='H')
                              a2[10][i]+=1;
                             else if( str.charAt(i)=='R')
                              a2[11][i]+=1;
                             else  if( str.charAt(i)=='N')
                            a2[12][i]+=1;
                             else if( str.charAt(i)=='Q')
                             a2[13][i]+=1;
                             else  if( str.charAt(i)=='E')
                             a2[14][i]+=1;
                             else if( str.charAt(i)=='D')
                             a2[15][i]+=1;
                             else  if( str.charAt(i)=='F')
                            a2[16][i]+=1;
                             else  if( str.charAt(i)=='W')
                            a2[17][i]+=1;
                             else  if( str.charAt(i)=='Y')
                             a2[18][i]+=1;
                             else if ( str.charAt(i)=='K')
                             a2[19][i]+=1;
                            
                             }
                             
                             }
           
       }
           
      
       }
        public static void freqcal(int array[][]) throws IOException    {
           
        for(int r=0;r<ROWS;r++)
        {
           
            for(int c=0;c<COLS;c++)
            {  
               double temp= array[r][c];
               
                temp=(temp/count);
              temp=roundtwofour( temp);
                af[r][c]=temp;
                
               
            }
        }
            
        
            
           
    }
        
         public static double scoring(String stri) 
    {
      
            double fscore=0.0;
               for(int i=0;i<stri.length();i++){
                  
                        if( stri.charAt(i)=='A')
                            fscore+=af[0][i];
                            else   if( stri.charAt(i)=='G')
                              fscore+=af[1][i];
                            else  if( stri.charAt(i)=='V')
                              fscore+=af[2][i];
                            else if( stri.charAt(i)=='L')
                              fscore+=af[3][i];
                            else if( stri.charAt(i)=='I')
                             fscore+=af[4][i];
                             else if( stri.charAt(i)=='S')
                             fscore+=af[5][i];
                             else if( stri.charAt(i)=='T')
                             fscore+=af[6][i];
                             else if( stri.charAt(i)=='C')
                            fscore+=af[7][i];
                             else if( stri.charAt(i)=='M')
                            fscore+=af[8][i];
                             else  if( stri.charAt(i)=='P')
                             fscore+=af[9][i];
                             else   if( stri.charAt(i)=='H')
                             fscore+=af[10][i];
                             else if( stri.charAt(i)=='R')
                              fscore+=af[11][i];
                             else  if( stri.charAt(i)=='N')
                            fscore+=af[12][i];
                             else if( stri.charAt(i)=='Q')
                             fscore+=af[13][i];
                             else  if( stri.charAt(i)=='E')
                             fscore+=af[14][i];
                             else if( stri.charAt(i)=='D')
                            fscore+=af[15][i];
                             else  if( stri.charAt(i)=='F')
                            fscore+=af[16][i];
                             else  if( stri.charAt(i)=='W')
                            fscore+=af[17][i];
                             else  if( stri.charAt(i)=='Y')
                             fscore+=af[18][i];
                             else if ( stri.charAt(i)=='K')
                             fscore+=af[19][i];
          }
           double cops=17.66;
         
           fscore=(fscore*cops);
         
           fscore=roundtwofour(fscore);
           return fscore;
   
       }
              
        
     public static double roundtwofour(double score)
    {
       // double d=Double.parseDouble(score);
        long l=(int) Math.round(score*10000);
        score=l/10000.0;
        return score;
    }  
      public static double  DimerScore (String stri) {
double dscore=0.0;        
double AA=0.4615;
double  AC=0.3684;
double AD=0.52;
double AE=0.5579;
double AF=0.551;
double AG=0.425;
double AH=0.5172;
double AI=0.54;
double AK=0.6309;
double AL=0.4706;
double AM=0.5417;
double AN=0.5833;
double AP=0.4118;
double AQ=0.4915;
double AR=0.519;
double AS=0.4719;
double AT=0.5;
double AV=0.5;
double AW=0.3846;
double AY=0.6471;
double CA = 0.3333 ;
double CC = 0.2308 ;
double CD = 0.3333 ;
double CE = 0.4 ;
double CF = 0.3077 ;
double CG = 0.2632 ;
double CH = 0.25 ;
double CI = 0.3077 ;
double CK = 0.375 ;
double CL = 0.2069 ;
double CM = 0.4 ;
double CN = 0.5 ;
double CP = 0.2105 ;
double CQ = 0.3125 ;
double CR = 0.4167 ;
double CS = 0.25 ;
double CT =0.3333 ;
double CV = 0.25 ;
double CY = 0.3333 ;
double DA = 0.569 ;
double DC = 0.3571 ;
double DD = 0.4615 ;
double DE = 0.5088 ;
double DF = 0.5476 ;
double DG = 0.5091 ;
double DH = 0.5 ;
double DI = 0.6034 ;
double DK = 0.6535 ;
double DL = 0.5319 ;
double DM = 0.4375 ;
double DN = 0.5714 ;
double DP = 0.5185 ;
double DQ = 0.5 ;
double DR = 0.5745 ;
double DS = 0.4237 ;
double DT = 0.4615 ;
double DV = 0.5192 ;
double DW = 0.5 ;
double DY = 0.6429 ;
double	EA=0.5392;
double	EC=0.3333;
double	ED=0.4932;
double	EE=0.5541;
double	EF=0.5909;
double	EG=0.5;
double	EH=0.6;
double	EI=0.5758;
double	EK=0.6386;
double	EL=0.5246;
double	EM=0.5667;
double	EN=0.5472;
double	EP=0.4151;
double	EQ=0.5179;
double	ER=0.5714;
double	ES=0.4783;
double	ET=0.4151;
double	EV=0.5455;
double	EW=0.4286;
double	EY=0.5882;
double	FA=0.6136;
double	FC=0.2308;
double	FD=0.6579;
double	FE=0.5952;
double	FF=0.5;
double	FG=0.6596;
double	FH=0.4762;
double	FI=0.5676;
double	FK=0.6267;
double	FL=0.4935;
double	FM=0.6667;
double	FN=0.5926;
double	FP=0.3333;
double	FQ=0.4706;
double	FR=0.55	;
double	FS=0.5246;
double	FT=0.4857;
double	FV=0.55	;
double	FW=0.4	;
double	FY=0.5417;
double	GA=0.4568;
double	GC=0.35	;
double	GD=0.46	;
double	GE=0.5068;
double	GF=0.5	;
double	GG=0.4494;
double	GH=0.4333;
double	GI=0.5385;
double	GK=0.6386;
double	GL=0.4545;
double	GM=0.5	;
double	GN=0.4884;
double	GP=0.3939;
double	GQ=0.4681;
double	GR=0.4595;
double	GS=0.4167;
double	GT=0.4821;
double	GV=0.5246;
double	GW=0.3571;
double	GY=0.5758;
double	HA=0.5769;
double	HC=0.2727;
double	HD=0.5625;
double	HE=0.5862;
double	HF=0.4545;
double	HG=0.5556;
double	HH=0.3077;
double	HI=0.4762;
double	HK=0.5581;
double	HL=0.5439;
double	HM=0.4444;
double	HN=0.5	;
double	HP=0.4643;
double	HQ=0.44	;
double	HR=0.5294;
double	HS=0.4324;
double	HT=0.3846;
double	HV=0.48	;
double	HW=0.4286;
double	HY=0.5333;
double	IA=0.6	;
double	IC=0.3333;
double	ID=0.5789;
double	IE=0.6441;
double	IF=0.5135;
double	IG=0.5897;
double	IH=0.5385;
double	II=0.5435;
double	IK=0.6538;
double	IL=0.5181;
double	IM=0.5882;
double	IN=0.5789;
double	IP=0.5	;
double	IQ=0.5122;
double	IR=0.6	;
double	IS=0.4677;
double	IT=0.5652;
double	IV=0.5882;
double	IW=0.4444;
double	IY=0.5517;
double	KA=0.5789;
double	KC=0.3333;
double	KD=0.6283;
double	KE=0.6023;
double	KF=0.7222;
double	KG=0.6015;
double	KH=0.661;
double	KI=0.6579;
double	KK=0.6429;
double	KL=0.598;
double	KM=0.6122;
double	KN=0.6444;
double	KP=0.5781;
double	KQ=0.5663;
double	KR=0.5167;
double	KS=0.5682;
double	KT=0.5727;
double	KV=0.5952;
double	KW=0.5455;
double	KY=0.7579;
double	LA=0.4732;
double	LC=0.2759;
double	LD=0.5	;
double	LE=0.5306;
double	LF=0.5085;
double	LG=0.4343;
double	LH=0.3913;
double	LI=0.493;
double	LK=0.6027;
double	LL=0.411;
double	LM=0.5	;
double	LN=0.5857;
double	LP=0.3605;
double	LQ=0.4565;
double	LR=0.5169;
double	LS=0.4496;
double	LT=0.4691;
double	LV=0.5	;
double	LW=0.3333;
double	LY=0.54	;
double	MA=0.44	;
double	MC=0.2	;
double	MD=0.619;
double	ME=0.5357;
double	MF=0.4545;
double	MG=0.5	;
double	MH=0.4444;
double	MI=0.5385;
double	MK=0.617;
double	ML=0.4516;
double	MM=0.625;
double	MN=0.4167;
double	MP=0.5	;
double	MQ=0.4667;
double	MR=0.5652;
double	MS=0.5	;
double	MT=0.5625;
double	MV=0.5217;
double	MY=0.5	;
double	NA=0.575;
double	NC=0.3636;
double	ND=0.5714;
double	NE=0.6087;
double	NF=0.5938;
double	NG=0.5476;
double	NH=0.5294;
double	NI=0.525;
double	NK=0.6364;
double	NL=0.5429;
double	NM=0.5385;
double	NN=0.5714;
double	NP=0.5263;
double	NQ=0.5	;
double	NR=0.5294;
double	NS=0.4681;
double	NT=0.5588;
double	NV=0.5952;
double	NW=0.5	;
double	NY=0.5238;
double	PA=0.4103;
double	PC=0.2143;
double	PD=0.4681;
double	PE=0.4861;
double	PF=0.5641;
double	PG=0.3837;
double	PH=0.4167;
double	PI=0.5455;
double	PK=0.5	;
double	PL=0.4	;
double	PM=0.4667;
double	PN=0.5172;
double	PP=0.3253;
double	PQ=0.3778;
double	PR=0.3934;
double	PS=0.3625;
double	PT=0.4043;
double	PV=0.5	;
double	PW=0.3571;
double	PY=0.5385;
double	QA=0.4444;
double	QC=0.25	;
double	QD=0.4474;
double	QE=0.4928;
double	QF=0.5	;
double	QG=0.4375;
double	QH=0.4	;
double	QI=0.6098;
double	QK=0.5688;
double	QL=0.475;
double	QM=0.4706;
double	QN=0.3929;
double	QP=0.3659;
double	QQ=0.451;
double	QR=0.5	;
double	QS=0.4186;
double	QT=0.4474;
double	QV=0.4783;
double	QW=0.4	;
double	QY=0.5909;
double	RA=0.5125;
double	RC=0.45	;
double	RD=0.5556;
double	RE=0.4872;
double	RF=0.561;
double	RG=0.4265;
double	RH=0.5806;
double	RI=0.6182;
double	RK=0.4628;
double	RL=0.5299;
double	RM=0.5652;
double	RN=0.5918;
double	RP=0.32	;
double	RQ=0.5273;
double	RR=0.4362;
double	RS=0.4875;
double	RT=0.5593;
double	RV=0.5652;
double	RW=0.4667;
double	RY=0.5625;
double	SA=0.475;
double	SC=0.2	;
double	SD=0.4909;
double	SE=0.4756;
double	SF=0.413;
double	SG=0.4725;
double	SH=0.3448;
double	SI=0.5	;
double	SK=0.5933;
double	SL=0.4331;
double	SM=0.4783;
double	SN=0.4889;
double	SP=0.3929;
double	SQ=0.4237;
double	SR=0.4079;
double	SS=0.3968;
double	ST=0.4762;
double	SV=0.4595;
double	SW=0.3158;
double	SY=0.5	;
double	TA=0.4833;
double	TC=0.2941;
double	TD=0.5122;
double	TE=0.5167;
double	TF=0.5	;
double	TG=0.5217;
double	TH=0.4286;
double	TI=0.5556;
double	TK=0.59	;
double	TL=0.4767;
double	TM=0.4444;
double	TN=0.5	;
double	TP=0.4576;
double	TQ=0.4211;
double	TR=0.5306;
double	TS=0.4714;
double	TT=0.4565;
double	TV=0.4407;
double	TW=0.3846;
double	TY=0.5	;
double	VA=0.4932;
double	VC=0.3182;
double	VD=0.6071;
double	VE=0.5769;
double	VF=0.5349;
double	VG=0.4746;
double	VH=0.5333;
double	VI=0.4898;
double	VK=0.6385;
double	VL=0.4857;
double	VM=0.4118;
double	VN=0.5581;
double	VP=0.5	;
double	VQ=0.4884;
double	VR=0.5574;
double	VS=0.411;
double	VT=0.5082;
double	VV=0.5072;
double	VW=0.375;
double	VY=0.5185;
double	WA=0.4	;
double	WC=0.25	;
double	WD=0.4545;
double	WE=0.4	;
double	WF=0.4286;
double	WG=0.3333;
double	WH=0.6667;
double	WI=0.3333;
double	WK=0.4643;
double	WL=0.4545;
double	WN=0.375;
double	WP=0.3333;
double	WQ=0.4	;
double	WR=0.3846;
double	WS=0.2727;
double	WT=0.4167;
double	WV=0.4167;
double	WW=0.25	;
double	WY=0.3333;
double	YA=0.6667;
double	YC=0.5	;
double	YD=0.6563;
double	YE=0.5758;
double	YF=0.5862;
double	YG=0.6667;
double	YH=0.5714;
double	YI=0.625;
double	YK=0.6769;
double	YL=0.4565;
double	YM=0.5556;
double	YN=0.64	;
double	YP=0.4286;
double	YQ=0.5	;
double	YR=0.5625;
double	YS=0.4865;
double	YT=0.6129;
double	YV=0.5938;
double	YW=0.4	;
double	YY=0.619;

       
  for(int i=0;i<stri.length()-1;i++)
        {
            if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='A')
            dscore+=AA; 
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='C')
            dscore+=AC;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='D')
            dscore+=AD;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='E')
            dscore+=AE;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='F')
            dscore+=AF;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='G')
            dscore+=AG;
           else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='H')
            dscore+=AH;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='I')
            dscore+=AI;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='K')
            dscore+=AK;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='L')
            dscore+=AL;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='M')
            dscore+=AM;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='N')
            dscore+=AN;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='P')
            dscore+=AP;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='Q')
            dscore+=AQ;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='R')
            dscore+=AR;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='S')
            dscore+=AS;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='T')
            dscore+=AT;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='V')
            dscore+=AV;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='W')
            dscore+=AW;
            else if(stri.charAt(i)=='A'&&stri.charAt(i+1)=='Y')
            dscore+=AY;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='A')
            dscore+=CA; 
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='C')
            dscore+=CC;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='D')
            dscore+=CD;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='E')
            dscore+=CE;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='F')
            dscore+=CF;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='G')
            dscore+=CG;
           else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='H')
            dscore+=CH;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='I')
            dscore+=CI;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='K')
            dscore+=CK;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='L')
            dscore+=CL;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='M')
            dscore+=CM;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='N')
            dscore+=CN;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='P')
            dscore+=CP;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='Q')
            dscore+=CQ;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='R')
            dscore+=CR;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='S')
            dscore+=CS;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='T')
            dscore+=CT;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='V')
            dscore+=CV;
            else if(stri.charAt(i)=='C'&&stri.charAt(i+1)=='Y')
            dscore+=CY;
            
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='A')
            dscore+=DA; 
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='C')
            dscore+=DC;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='D')
            dscore+=DD;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='E')
            dscore+=DE;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='F')
            dscore+=DF;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='G')
            dscore+=DG;
           else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='H')
            dscore+=DH;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='I')
            dscore+=DI;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='K')
            dscore+=DK;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='L')
            dscore+=DL;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='M')
            dscore+=DM;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='N')
            dscore+=DN;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='P')
            dscore+=DP;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='Q')
            dscore+=DQ;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='R')
            dscore+=DR;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='S')
            dscore+=DS;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='T')
            dscore+=DT;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='V')
            dscore+=DV;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='W')
            dscore+=DW;
            else if(stri.charAt(i)=='D'&&stri.charAt(i+1)=='Y')
            dscore+=DY;
            
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='A')
            dscore+=EA; 
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='C')
            dscore+=EC;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='D')
            dscore+=ED;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='E')
            dscore+=EE;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='F')
            dscore+=EF;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='G')
            dscore+=EG;
           else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='H')
            dscore+=EH;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='I')
            dscore+=EI;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='K')
            dscore+=EK;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='L')
            dscore+=EL;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='M')
            dscore+=EM;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='N')
            dscore+=EN;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='P')
            dscore+=EP;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='Q')
            dscore+=EQ;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='R')
            dscore+=ER;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='S')
            dscore+=ES;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='T')
            dscore+=ET;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='V')
            dscore+=EV;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='W')
            dscore+=EW;
            else if(stri.charAt(i)=='E'&&stri.charAt(i+1)=='Y')
            dscore+=EY;
            
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='A')
            dscore+=FA; 
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='C')
            dscore+=FC;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='D')
            dscore+=FD;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='E')
            dscore+=FE;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='F')
            dscore+=FF;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='G')
            dscore+=FG;
           else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='H')
            dscore+=FH;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='I')
            dscore+=FI;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='K')
            dscore+=FK;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='L')
            dscore+=FL;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='M')
            dscore+=FM;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='N')
            dscore+=FN;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='P')
            dscore+=FP;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='Q')
            dscore+=FQ;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='R')
            dscore+=FR;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='S')
            dscore+=FS;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='T')
            dscore+=FT;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='V')
            dscore+=FV;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='W')
            dscore+=FW;
            else if(stri.charAt(i)=='F'&&stri.charAt(i+1)=='Y')
            dscore+=FY;
            
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='A')
            dscore+=GA; 
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='C')
            dscore+=GC;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='D')
            dscore+=GD;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='E')
            dscore+=GE;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='F')
            dscore+=GF;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='G')
            dscore+=GG;
           else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='H')
            dscore+=GH;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='I')
            dscore+=GI;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='K')
            dscore+=GK;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='L')
            dscore+=GL;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='M')
            dscore+=GM;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='N')
            dscore+=GN;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='P')
            dscore+=GP;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='Q')
            dscore+=GQ;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='R')
            dscore+=GR;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='S')
            dscore+=GS;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='T')
            dscore+=GT;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='V')
            dscore+=GV;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='W')
            dscore+=GW;
            else if(stri.charAt(i)=='G'&&stri.charAt(i+1)=='Y')
            dscore+=GY;
                        
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='A')
            dscore+=HA; 
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='C')
            dscore+=HC;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='D')
            dscore+=HD;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='E')
            dscore+=HE;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='F')
            dscore+=HF;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='G')
            dscore+=HG;
           else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='H')
            dscore+=HH;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='I')
            dscore+=HI;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='K')
            dscore+=HK;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='L')
            dscore+=HL;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='M')
            dscore+=HM;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='N')
            dscore+=HN;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='P')
            dscore+=HP;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='Q')
            dscore+=HQ;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='R')
            dscore+=HR;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='S')
            dscore+=HS;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='T')
            dscore+=HT;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='V')
            dscore+=HV;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='W')
            dscore+=HW;
            else if(stri.charAt(i)=='H'&&stri.charAt(i+1)=='Y')
            dscore+=HY;
            
            
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='A')
            dscore+=IA; 
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='C')
            dscore+=IC;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='D')
            dscore+=ID;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='E')
            dscore+=IE;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='F')
            dscore+=IF;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='G')
            dscore+=IG;
           else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='H')
            dscore+=IH;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='I')
            dscore+=II;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='K')
            dscore+=IK;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='L')
            dscore+=IL;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='M')
            dscore+=IM;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='N')
            dscore+=IN;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='P')
            dscore+=IP;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='Q')
            dscore+=IQ;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='R')
            dscore+=IR;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='S')
            dscore+=IS;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='T')
            dscore+=IT;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='V')
            dscore+=IV;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='W')
            dscore+=IW;
            else if(stri.charAt(i)=='I'&&stri.charAt(i+1)=='Y')
            dscore+=IY;
           
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='A')
            dscore+=KA; 
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='C')
            dscore+=KC;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='D')
            dscore+=KD;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='E')
            dscore+=KE;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='F')
            dscore+=KF;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='G')
            dscore+=KG;
           else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='H')
            dscore+=KH;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='I')
            dscore+=KI;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='K')
            dscore+=KK;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='L')
            dscore+=KL;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='M')
            dscore+=KM;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='N')
            dscore+=KN;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='P')
            dscore+=KP;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='Q')
            dscore+=KQ;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='R')
            dscore+=KR;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='S')
            dscore+=KS;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='T')
            dscore+=KT;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='V')
            dscore+=KV;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='W')
            dscore+=KW;
            else if(stri.charAt(i)=='K'&&stri.charAt(i+1)=='Y')
            dscore+=KY;
   
            
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='A')
            dscore+=LA; 
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='C')
            dscore+=LC;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='D')
            dscore+=LD;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='E')
            dscore+=LE;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='F')
            dscore+=LF;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='G')
            dscore+=LG;
           else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='H')
            dscore+=LH;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='I')
            dscore+=LI;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='K')
            dscore+=LK;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='L')
            dscore+=LL;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='M')
            dscore+=LM;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='N')
            dscore+=LN;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='P')
            dscore+=LP;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='Q')
            dscore+=LQ;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='R')
            dscore+=LR;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='S')
            dscore+=LS;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='T')
            dscore+=LT;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='V')
            dscore+=LV;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='W')
            dscore+=LW;
            else if(stri.charAt(i)=='L'&&stri.charAt(i+1)=='Y')
            dscore+=LY;
               
              else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='A')
            dscore+=MA; 
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='C')
            dscore+=MC;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='D')
            dscore+=MD;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='E')
            dscore+=ME;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='F')
            dscore+=MF;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='G')
            dscore+=MG;
           else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='H')
            dscore+=MH;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='I')
            dscore+=MI;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='K')
            dscore+=MK;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='L')
            dscore+=ML;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='M')
            dscore+=MM;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='N')
            dscore+=MN;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='P')
            dscore+=MP;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='Q')
            dscore+=MQ;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='R')
            dscore+=MR;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='S')
            dscore+=MS;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='T')
            dscore+=MT;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='V')
            dscore+=MV;
            else if(stri.charAt(i)=='M'&&stri.charAt(i+1)=='Y')
            dscore+=MY;
            
                        
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='A')
            dscore+=NA; 
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='C')
            dscore+=NC;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='D')
            dscore+=ND;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='E')
            dscore+=NE;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='F')
            dscore+=NF;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='G')
            dscore+=NG;
           else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='H')
            dscore+=NH;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='I')
            dscore+=NI;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='K')
            dscore+=NK;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='L')
            dscore+=NL;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='M')
            dscore+=NM;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='N')
            dscore+=NN;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='P')
            dscore+=NP;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='Q')
            dscore+=NQ;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='R')
            dscore+=NR;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='S')
            dscore+=NS;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='T')
            dscore+=NT;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='V')
            dscore+=NV;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='W')
            dscore+=NW;
            else if(stri.charAt(i)=='N'&&stri.charAt(i+1)=='Y')
            dscore+=NY;
            
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='A')
            dscore+=PA; 
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='C')
            dscore+=PC;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='D')
            dscore+=PD;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='E')
            dscore+=PE;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='F')
            dscore+=PF;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='G')
            dscore+=PG;
           else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='H')
            dscore+=PH;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='I')
            dscore+=PI;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='K')
            dscore+=PK;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='L')
            dscore+=PL;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='M')
            dscore+=PM;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='N')
            dscore+=PN;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='P')
            dscore+=PP;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='Q')
            dscore+=PQ;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='R')
            dscore+=PR;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='S')
            dscore+=PS;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='T')
            dscore+=PT;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='V')
            dscore+=PV;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='W')
            dscore+=PW;
            else if(stri.charAt(i)=='P'&&stri.charAt(i+1)=='Y')
            dscore+=PY;
            
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='A')
            dscore+=QA; 
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='C')
            dscore+=QC;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='D')
            dscore+=QD;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='E')
            dscore+=QE;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='F')
            dscore+=QF;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='G')
            dscore+=QG;
           else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='H')
            dscore+=QH;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='I')
            dscore+=QI;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='K')
            dscore+=QK;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='L')
            dscore+=QL;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='M')
            dscore+=QM;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='N')
            dscore+=QN;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='P')
            dscore+=QP;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='Q')
            dscore+=QQ;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='R')
            dscore+=QR;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='S')
            dscore+=QS;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='T')
            dscore+=QT;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='V')
            dscore+=QV;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='W')
            dscore+=QW;
            else if(stri.charAt(i)=='Q'&&stri.charAt(i+1)=='Y')
            dscore+=QY;
               
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='A')
            dscore+=RA; 
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='C')
            dscore+=RC;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='D')
            dscore+=RD;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='E')
            dscore+=RE;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='F')
            dscore+=RF;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='G')
            dscore+=RG;
           else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='H')
            dscore+=RH;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='I')
            dscore+=RI;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='K')
            dscore+=RK;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='L')
            dscore+=RL;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='M')
            dscore+=RM;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='N')
            dscore+=RN;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='P')
            dscore+=RP;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='Q')
            dscore+=RQ;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='R')
            dscore+=RR;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='S')
            dscore+=RS;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='T')
            dscore+=RT;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='V')
            dscore+=RV;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='W')
            dscore+=RW;
            else if(stri.charAt(i)=='R'&&stri.charAt(i+1)=='Y')
            dscore+=RY;
            
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='A')
            dscore+=SA; 
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='C')
            dscore+=SC;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='D')
            dscore+=SD;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='E')
            dscore+=SE;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='F')
            dscore+=SF;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='G')
            dscore+=SG;
           else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='H')
            dscore+=SH;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='I')
            dscore+=SI;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='K')
            dscore+=SK;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='L')
            dscore+=SL;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='M')
            dscore+=SM;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='N')
            dscore+=SN;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='P')
            dscore+=SP;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='Q')
            dscore+=SQ;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='R')
            dscore+=SR;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='S')
            dscore+=SS;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='T')
            dscore+=ST;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='V')
            dscore+=SV;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='W')
            dscore+=SW;
            else if(stri.charAt(i)=='S'&&stri.charAt(i+1)=='Y')
            dscore+=SY;
           
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='A')
            dscore+=TA; 
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='C')
            dscore+=TC;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='D')
            dscore+=TD;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='E')
            dscore+=TE;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='F')
            dscore+=TF;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='G')
            dscore+=TG;
           else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='H')
            dscore+=TH;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='I')
            dscore+=TI;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='K')
            dscore+=TK;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='L')
            dscore+=TL;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='M')
            dscore+=TM;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='N')
            dscore+=TN;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='P')
            dscore+=TP;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='Q')
            dscore+=TQ;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='R')
            dscore+=TR;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='S')
            dscore+=TS;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='T')
            dscore+=TT;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='V')
            dscore+=TV;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='W')
            dscore+=TW;
            else if(stri.charAt(i)=='T'&&stri.charAt(i+1)=='Y')
            dscore+=TY;
            
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='A')
            dscore+=VA; 
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='C')
            dscore+=VC;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='D')
            dscore+=VD;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='E')
            dscore+=VE;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='F')
            dscore+=VF;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='G')
            dscore+=VG;
           else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='H')
            dscore+=VH;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='I')
            dscore+=VI;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='K')
            dscore+=VK;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='L')
            dscore+=VL;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='M')
            dscore+=VM;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='N')
            dscore+=VN;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='P')
            dscore+=VP;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='Q')
            dscore+=VQ;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='R')
            dscore+=VR;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='S')
            dscore+=VS;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='T')
             dscore+=VT;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='V')
            dscore+=VV;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='W')
            dscore+=VW;
            else if(stri.charAt(i)=='V'&&stri.charAt(i+1)=='Y')
            dscore+=VY;
            
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='A')
            dscore+=WA; 
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='C')
            dscore+=WC;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='D')
            dscore+=WD;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='E')
            dscore+=WE;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='F')
            dscore+=WF;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='G')
            dscore+=WG;
           else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='H')
            dscore+=WH;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='I')
            dscore+=WI;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='K')
            dscore+=WK;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='L')
            dscore+=WL;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='N')
            dscore+=WN;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='P')
            dscore+=WP;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='Q')
            dscore+=WQ;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='R')
            dscore+=WR;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='S')
            dscore+=WS;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='T')
            dscore+=WT;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='V')
            dscore+=WV;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='W')
            dscore+=WW;
            else if(stri.charAt(i)=='W'&&stri.charAt(i+1)=='Y')
            dscore+=WY;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='A')
            dscore+=YA; 
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='C')
            dscore+=YC;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='D')
            dscore+=YD;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='E')
            dscore+=YE;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='F')
            dscore+=YF;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='G')
            dscore+=YG;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='H')
            dscore+=YH;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='I')
            dscore+=YI;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='K')
            dscore+=YK;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='L')
            dscore+=YL;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='N')
            dscore+=YN;
             else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='M')
            dscore+=YM;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='P')
            dscore+=YP;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='Q')
            dscore+=YQ;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='R')
            dscore+=YR;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='S')
            dscore+=YS;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='T')
            dscore+=YT;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='V')
            dscore+=YV;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='W')
            dscore+=YW;
            else if(stri.charAt(i)=='Y'&&stri.charAt(i+1)=='Y')
            dscore+=YY;
            }
          double codi=5.553;
          dscore=(dscore*codi);
          return dscore;
      }



}







