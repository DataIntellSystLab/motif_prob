#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <bits/stdc++.h>
#include <cctype>
#include <ctime>
#include <vector>
#include <unistd.h>  
#include <ctype.h>

using namespace std;

//global variables
long double S0[500],S3[500],SR[500],COMB[5000],tabR[5000],tabF[5000];
string tabmsg[4];


//calculating ik-th element of S0 and S3
 int ROUTCALCS3(int IK, long long int k, long long int m, long long int n, long double PB, long double maxld)
 {
   long double limovflw1,limovflw2;
   int riscalcS3=0;
   limovflw1=maxld / k;
   limovflw2=maxld /S0[IK-m];

   if ((S0[IK-1] >= limovflw1) || (PB >= limovflw2))
      {
       tabmsg[0]="******************************************";
       tabmsg[1]="* ERROR IN THE CALCULATION OF S0         *";
       tabmsg[2]="* error type  NUMERICAL CAPACITY  index IK: " + to_string(IK);
       tabmsg[3]="*******      O W E R F L O W       *******";
       riscalcS3=99;
       return riscalcS3;
      }
   S0[IK]=(S0[IK-1] * k) - (PB * S0[IK-m]);
   if (isinf(S0[IK]))
      {
       tabmsg[0]="******************************************";
       tabmsg[1]="* ERROR IN THE CALCULATION OF S0         *";
       tabmsg[2]="* error type INF  index IK: " + to_string(IK);
       tabmsg[3]="******************************************";
       riscalcS3=99;
       return riscalcS3;
      }
   if (isnan(S0[IK]))
      {
       tabmsg[0]="******************************************";
       tabmsg[1]="* ERROR IN THE CALCULATION OF S0         *";
       tabmsg[2]="* error type NAN  index IK: " + to_string(IK);
       tabmsg[3]="******************************************";
       riscalcS3=99;
       return riscalcS3;
      }
   S3[IK]=S0[IK]/(pow(k,IK));
   if (isinf(S3[IK]))
      {
       tabmsg[0]="******************************************";
       tabmsg[1]="* ERROR IN THE CALCULATION OF S3         *";
       tabmsg[2]="* error type INF  index IK: " + to_string(IK);
       tabmsg[3]="******************************************";
       riscalcS3=99;
       return riscalcS3;
      }
   if (isnan(S3[IK]))
      {
       tabmsg[0]="******************************************";
       tabmsg[1]="* ERROR IN THE CALCULATION OF S3         *";
       tabmsg[2]="* error type NAN  index IK: " + to_string(IK);
       tabmsg[3]="******************************************";
       riscalcS3=99;
       return riscalcS3;
      }
    return riscalcS3;
  }


// calculating tabR ik-th element power
int ROUTCALCPOWR(int IK, long long int m, long long int n, long long int nepsilon, long double minSR, long double maxld)
{
    long double E=2.71828182845904525536;
    long double y,z,limovflw1;
    int riscalcpowr=0;
    y=(n - (m*IK) -  nepsilon * (IK + 1)) * log10(minSR);
    z=log(10);
    limovflw1=log(maxld);
    if ((y*z) >= limovflw1)
       {
         tabmsg[0]="******************************************";
         tabmsg[1]="* ERROR IN CALCULATION OF R POWER        *";
         tabmsg[2]="* error type NUMERICAL CAPACITY  index IK: " + to_string(IK);
         tabmsg[3]="*******      O W E R F L O W       *******";
         riscalcpowr=99;
         return riscalcpowr;
        }
    tabR[IK]=pow(E,(y*z));
    if (isinf(tabR[IK]))
       {
         tabmsg[0]="******************************************";
         tabmsg[1]="* ERROR IN CALCULATION OF R POWER        *";
         tabmsg[2]="* error type INF  index IK: " + to_string(IK);
         tabmsg[3]="******************************************";
         riscalcpowr=99;
         return riscalcpowr;
       }
     if (isnan(tabR[IK]))
       {
         tabmsg[0]="******************************************";
         tabmsg[1]="* ERROR IN CALCULATION OF R POWER        *";
         tabmsg[2]="* error type NAN  index IK: " + to_string(IK);
         tabmsg[3]="******************************************";
         riscalcpowr=99;
         return riscalcpowr;
       }
    return riscalcpowr;
}


//calcularing COMB ik-th element
int ROUTCALCCOMB(int IK, long long int m, long long int n, long double maxld)
{
  long double  numeratore,denominatore;
  long long int xn;
  long double limovflw1;
  int riscalccomb=0;
  xn=n+IK-m*IK;
  numeratore=xn;
  denominatore=IK;
  for(int i=1;i<IK;i++)
     {
       limovflw1=maxld /(xn-i);
       if (numeratore >= limovflw1)
          {
           tabmsg[0]="********************************************";
           tabmsg[1]="* ERROR IN CALCULATION OF THE COMBINATIONS *";
           tabmsg[2]="* error type NUMERICAL CAPACITY  index IK: " + to_string(IK);
           tabmsg[3]="*******      O W E R F L O W       *******";
           riscalccomb=99;
           return riscalccomb;
          }
        numeratore=numeratore*(xn-i);
        if (isinf(numeratore))
           {
            tabmsg[0]="********************************************";
            tabmsg[1]="* ERROR IN CALCULATION OF THE COMBINATIONS *";
            tabmsg[2]="* error type INF index IK: " + to_string(IK);
            tabmsg[3]="********************************************";
            riscalccomb=99;
            return riscalccomb;
           }
        if (isnan(numeratore))
           {
             tabmsg[0]="********************************************";
             tabmsg[1]="* ERROR IN CALCULATION OF THE COMBINATIONS *";
             tabmsg[2]="* error type NAN index IK: " + to_string(IK);
             tabmsg[3]="********************************************";
             riscalccomb=99;
             return riscalccomb;
           }
         denominatore=denominatore*(IK-i);
       }
       COMB[IK]=numeratore/denominatore;
       return riscalccomb;
  }

 static void show_usage(std::string name)
 {
   std::cerr << "\nUsage: " << name << " <option(s)>\n"
     << "Options:\n"
     << " -i INPUT_FILE -o OUTPUT_FILE {-f GENOME_FILE | -a freq_a -c freq_c -g freq_g -t freq_t -n genome_length}\n"
     << "-i input file\n"
     << "-o output file\n"
     << "-a A frequency\n"
     << "-c C frequency\n"
     << "-g G frequency\n"
     << "-t T frequency\n"
     << "-n total genome length\n"
     << "-f genome file. If specfied, nucleotide frequencies will be calculated on this file, and other options will be ignored.\n"
     << "   Fasta format. Note that the total genome length will be calculated on the fasta full sequence lenght, while nucleotide\n"
     << "   frequencies will only consider (A, C, G, T) characters.\n"
     << "Examples:\n formlp03 -i INPCPLUSPLUS_unix.TXT -o OUTCPLUSPLUS.TXT -f mygenom.fasta\n formlp03 -i INPCPLUSPLUS_unix.TXT -o OUTCPLUSPLUS.TXT -a 0.1 -c 0.3 -g 0.1 -t 0.5 -n 100000"
     << std::endl;
 }
 

int main(int argc, char* argv[])
{
  int c;
  char *infile = NULL;
  char *outfile = NULL;
  char *genomefile = NULL;
  long long int k,m,n,l,indminSR,nepsilon,jlim;  
  long double pA,pC,pG,pT,pTOT,PS,PB,eps;
  n = pA = pC = pG = pT = -1; // if user won't provide values, will be calculated from genome file
  eps = 0.0001; // tolerance for probability sum check
  
  
  while ((c = getopt(argc, argv, "i:o:a:c:g:t:n:f:")) != -1)
    switch (c) {
    case 'i':
      infile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'a':
      pA = stod(optarg);
      break;
    case 'c':
      pC = stod(optarg);
      break;
    case 'g':
      pG = stod(optarg);
      break;
    case 't':
      pT = stod(optarg);
      break;
    case 'n':
      n = stoll(optarg);
      break;
    case 'f':
      genomefile = optarg;
      break;
    case '?':
      printf("\nWrong option configuration.");
      show_usage(argv[0]);
      return 1;
    case ':':
      printf("\nWrong option configuration.");
      show_usage(argv[0]);
      return 1;
    }
  
 
  pTOT=pA+pC+pG+pT;

  if(genomefile != NULL){	// option -f present, calculate frequencies from provided genome fasta file
    if(access( genomefile, F_OK ) != 0 ) {
      std::cerr << "\nGenome file not found. Genome file searched: "
                << genomefile
                << std::endl;
      show_usage(argv[0]);
      return 1;
    }
    else{
      string nomegenomefile = genomefile;
      fstream ginp(nomegenomefile,fstream::in);
    
      std::string seq;
      unordered_map<char, unsigned long long int> Freqs; 
      n=0; 
      while(getline(ginp, seq)){
        if(seq.at(0) != '>'){
          std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
          n+=seq.length(); 
            for (int i = 0; seq[i]; i++){
            if (Freqs.find(seq[i]) == Freqs.end()){
              Freqs.insert(make_pair(seq[i], 1));
            }
            else{
              Freqs[seq[i]]++;          
            }
          }
        }
      }
      std::cerr << "Counted in genome file " << genomefile << std::endl;        
      for (auto& it : Freqs) {
        std::cerr << it.first << " = " << it.second << std::endl;  
      }
      long long int total_standard_nucleotides = Freqs['A']+Freqs['C']+Freqs['G']+Freqs['T'];
      pA=Freqs['A']/float(total_standard_nucleotides);
      pC=Freqs['C']/float(total_standard_nucleotides);
      pG=Freqs['G']/float(total_standard_nucleotides);
      pT=Freqs['T']/float(total_standard_nucleotides);
    }
  }
  else{
    if(pTOT < 1-eps or pTOT > 1+eps or pA < 0 or pC < 0 or pG < 0 or pT < 0){
      std::cerr << "\nProblem with the options. Sum of the nucleotide frequencies should be 1.\n"
                << std::endl;
      show_usage(argv[0]);
      return 1;
    }
    
    if(pA < 0 or pC < 0 or pG < 0 or pT < 0 or n < 0){
      std::cerr << "\nProblem with the options. Probabilities and genome length should positive numbers.\n"
                << std::endl;
      show_usage(argv[0]);
      return 1;
    }

  }

  string nomefilein="INPCPLUSPLUS.TXT";
  if(infile == NULL){
    std::cerr << "\nInput file not specificed, proceeding with default one:\n"
              << nomefilein
              << std::endl;
  }else{
    if(access( infile, F_OK ) != 0 ) {
      std::cerr << "\nInput file not found."
                << std::endl;
      show_usage(argv[0]);
      return 1;
    }
    
    nomefilein.assign(infile);
    }

  string nomefileout="OUTCPLUSPLUS.TXT";
  if(outfile == NULL) {
    std::cerr << "\nOutput file not specificed, proceeding with default one:\n"
              << nomefileout
              << std::endl;
  } else{
    nomefileout.assign(outfile);
    }

  fstream stat(nomefileout,fstream::out);

  long double  minSR,PS0NE,SOMF,PAR1,PAR2,PAR3,PAR4;
  float fpa,fpc,fpg,fpt,fptot;

  long double limite100;
  int tappoparz, letture=0, stringheesatte=0, stringheanomale=0;
  long double  maxld;
  maxld=numeric_limits<long double>::max();
  string xrecord, stringa;

  fstream finp(nomefilein,fstream::in);
  k=4;

  //INIZIO
  time_t start = time(0), stop;
  tappoparz=500;
  
  stat<<"   CALCULATION OF PROBABILITY DISTRIBUTION "<<endl;
  stat<<"number of characters ............ "<<k<<endl;
  stat<<"sequence length ................. "<<n<<endl;
  stat<<"Probability A ................... "<<pA<<endl;
  stat<<"Probability C ................... "<<pC<<endl;
  stat<<"Probability G ................... "<<pG<<endl;
  stat<<"Probability T ................... "<<pT<<endl;
  stat<<"number of cycles of J ........... "<<tappoparz<<endl;
  stat<<"Output file name ................ "<<nomefileout<<endl;
  stat<<"Input file name ................. "<<nomefilein<<endl;
  stat<<endl;
  
  if (finp)
     {
       while(getline(finp,xrecord))
       {   stringa=xrecord;
  
  // CALCOLAPROBABILITASTRINGA:
       int swsalta=0;
       PS=1;
       m=0;
       l=stringa.size();
       stat<<"String              ............. "<<stringa<<endl;
       for (int i=0;i<l;i++)
       {
         string elem;
         elem=stringa.substr(i,1);
         if (elem=="A"  or elem=="C" or elem=="G" or elem=="T")
            {
              if (elem == "A")
                 PS=PS*pA;
              if (elem == "C")
                 PS=PS*pC;
              if (elem == "G")
                 PS=PS*pG;
              if (elem == "T")
                  PS=PS*pT;
               m++;
            }
         else
            {
              stat<<"**********************"<<endl;
              stat<<"* WRONG TYPED STRING *"<<endl;
              stat<<"**********************"<<endl;
              stat<<endl;
              i=l;
              swsalta=1;
            }
  
        }
      if (swsalta==1)
          continue;
      stat<<"string length   ................. "<<m<<endl;
      // CONTROLLOCLUMPIZZABILE:
      string xini,xfin;
      int swclump,tappoclump;
      swclump=0;
      tappoclump=m / 2;
      for (int i=1;i<tappoclump+1;i++)
      {
         xini=stringa.substr(0,i);
         xfin=stringa.substr(m-i,i);
         if (xini==xfin)
         {
             swclump=1;
             i=tappoclump+1;
         }
      }
         if (swclump==1)
             stat<<"Clumped string .................. "<<"Y"<<endl;
         else
             stat<<"Clumped string .................. "<<"N"<<endl;
  
       // Check lenghts:
       if (m>n)
         {
           stat<<"***********************************"<<endl;
           stat<<"* STRING LENGTH > SEQUENCE LENGTH *"<<endl;
           stat<<"***********************************"<<endl;
           stat<<endl;
           swsalta=1;
         }
          if (swsalta==1)
             continue;
  
  
  
      // CONTROLLODATI:
      PB=PS * pow(k,m);
  
      // CALCOLO J LIMITE
       jlim=n/m;
       if (jlim>tappoparz)
           jlim=tappoparz;
  
      stat.precision(25);
  
      stat<<"String probability PS ........... "<<PS<<endl;
      stat<<"Probability PS*k**m ............. "<<PB<<endl;
      stat<<"j limit calculated .............. "<<jlim<<endl;
  
      // calcolo limite probabilit? 100
      //limite100 = (n-m+1) * log(1.0-PS);
  
      //if (limite100 > -0.0001)
      //    stat<<"*** EXPECTED PROBABILITY 100% *** "<<endl;
  
      //***************
      //CALCOLO MIN R
      //     NEPSILON
      //***************
  
      minSR=1;
      indminSR=1;
      SR[m-1]=1;
      for(int i=0;i<m;i++)
      {
        S0[i]=pow(k,i);
        S3[i]=1;
      }
      S0[m]=(S0[m-1] * k) - (PB * S0[0]);
      S3[m]=S0[m]/(pow(k,m));
  
      for (int i=m;i<tappoparz;i++)
      {
          int IK=i+1;
          int riscalcS3;
          riscalcS3 = ROUTCALCS3(IK, k, m, n, PB, maxld);
          if (riscalcS3==99)
             {
              stat<<tabmsg[0]<<endl;
              stat<<tabmsg[1]<<endl;
              stat<<tabmsg[2]<<endl;
              stat<<tabmsg[3]<<endl;
              goto CONTINUA_LETTURA;
             }
  
          SR[i]=S3[i+1]/S3[i];
  
          if(SR[i]<minSR)
            {
              minSR=SR[i];
              indminSR=i;
            }
          if(SR[i]==SR[i-1])
               i=tappoparz;
      }
      if (indminSR==1)
         {
           nepsilon=1;
           PS0NE=S3[1];
         }
      else
         {
           nepsilon=indminSR-2;
           PS0NE=S3[nepsilon];
         }
  
  
      stat<<"minimum value of R .............. "<<minSR<<endl;
      stat<<"Index of minimum value .......... "<<indminSR<<endl;
      stat<<"n-epsilon value ................. "<<nepsilon<<endl;
      stat<<"PS(0,ne) value .................. "<<PS0NE<<endl;
      stat<<endl;
      //***************
      //***************
      //***************
  
  
  
  
  
      //***************
      // CALCOLO
      // FORMULA FINALE
      //***************
      double differenza;
      stat<<endl;
      stat<<"    FINAL FORMULA CALCULATION"<<endl;
      SOMF=0;
      COMB[0]=1;
      COMB[1]=n-m+1;
      int jstampa;
      for (int j=0;j<jlim;j++)
          {
            PAR1=(pow(PS,j));
            int IK=j;
            int riscalcpowr;
            riscalcpowr = ROUTCALCPOWR(IK, m, n, nepsilon, minSR, maxld);
            if (riscalcpowr==99)
               {
                stat<<tabmsg[0]<<endl;
                stat<<tabmsg[1]<<endl;
                stat<<tabmsg[2]<<endl;
                stat<<tabmsg[3]<<endl;
                goto CONTINUA_LETTURA;
               }
  
            PAR2=(tabR[j]);
            PAR3=(pow(PS0NE,(j+1)));
            if (j<2)
                goto SALTACALCCOMB;
            IK=j;
            int riscalccomb;
            riscalccomb = ROUTCALCCOMB(IK, m, n, maxld);
            if (riscalccomb==99)
               {
                stat<<tabmsg[0]<<endl;
                stat<<tabmsg[1]<<endl;
                stat<<tabmsg[2]<<endl;
                stat<<tabmsg[3]<<endl;
                goto CONTINUA_LETTURA;
               }
  
    SALTACALCCOMB:
            PAR4=(COMB[j]);
            tabF[j]=PAR1*PAR2*PAR3*PAR4*100;
            if (isinf(tabF[j]))
              {
                  stat<<"******************************************"<<endl;
                  stat<<"* ERROR IN CALCULATION OF THE FORMULA    *"<<endl;
                  stat<<"* error type INF  index j:" <<j<<endl;
                  stat<<"******************************************"<<endl;
                  goto CONTINUA_LETTURA;
              }
            if (isnan(tabF[j]))
              {
                  stat<<"******************************************"<<endl;
                  stat<<"* ERROR IN CALCULATION OF THE FORMULA    *"<<endl;
                  stat<<"* error type NAN  index j:" <<j<<endl;
                  stat<<"******************************************"<<endl;
                  goto CONTINUA_LETTURA;
              }
  
            SOMF=SOMF+tabF[j];
             if (tabF[j]>=0.000001)
                 stat<<"P("<<j<<","<<m<<","<<n<<") = "<<tabF[j]<<endl;
             else
                  {
                    differenza=100-SOMF;
                    if (differenza<0.000001)
                       {j=jlim;}
                  }
  
  
          }
       if (SOMF > 105)
           {stat<<"*** abnormal situation ***"<<endl;stringheanomale++;}
       if (SOMF < 95)
           {stat<<"*** abnormal situation ***"<<endl;stringheanomale++;}
       stringheesatte++;
  
      stat<<endl;
      //stat.precision(10);
      //stat<<"j medio        : "<<(PS * n)<<endl;
      stat<<"Sum probability: "<<SOMF<<endl;
  
      //***************
      //***************
      //***************
  
  CONTINUA_LETTURA:
       stat<<"============================"<<endl<<endl<<endl;
       }
     }
  
  stat<<"processed strings .... : "<<stringheesatte<<endl;
  stat<<"abnormal strings ...... : "<<stringheanomale<<endl;
  stat << endl;
  stop=time(0);
  stat << endl;
  stat << "start: " << start << endl;
  stat << "stop : " << stop << endl;
  int diff=stop - start;
  stat << "processing times sec : " << diff << endl;
  stat.close();
  finp.close();
  
  
  return 0;
}
