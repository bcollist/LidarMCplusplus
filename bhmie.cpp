#include "bhmie.hpp"
#include "constants.hpp"
#include "iostream"

using namespace std;

int bhmie(double x, complex<double> refrel, int nang, double* Qscat_p, double* Qext_p, double* Qabs_p, double* Qback_p, complex<double>* S1_p, complex<double>* S2_p){

    // Variable Definitions
    complex<double> y; double dx; double nstop; double ymod; int nmx; double dang; double theta;
    int nn;
    double RN; double DN; double FN;
    complex<double> AN; complex<double> BN;
    double PSI; double PSI0; double PSI1;
    double CHI; double CHI0; double CHI1;
    double APSI; double APSI0; double APSI1;
    complex <double> XI; complex <double> XI0; complex <double> XI1;
    double P; double T;

    // Assign variables that will be exported to the value of their pointers
    double Qscat_f = *Qscat_p; double Qext_f = *Qext_p; double Qabs_f = *Qabs_p; double Qback_f = *Qback_p;

    // Array Definitions
    double PI[nang+1]; double PI0[nang+1]; double PI1[nang+1];
    complex<double> S1[2*nang+1]; complex<double> S2[2*nang+1];
    double AMU[nang+1]; double TAU[nang+1];

    dx = x;
    y = x * refrel;

    nstop = ceil(x + 4 * pow(x,0.3333) +2);
    ymod = abs(y);
    nmx = max(nstop,ceil(ymod))+15;
    complex<double> D[nmx];
    dang = pi/2/(nang-1);


    for (int i = 1; i<=nang; i++){
        theta = (double)(i-1) * dang;
        AMU[i] =  cos(theta);
    }

    //Logarithmic derivative D(j) calculated by downward recurence
    D[nmx]= complex <double>(0,0);
    nn = nmx-1;

    for (int n = 1; n<=nn; n++){
        RN = nmx-n+1;
        D[nmx-n]=(RN/y)-(1.0/(D[nmx-n+1]+RN/y));
    }

    for (int j = 1; j <= nang; j++){
        PI0[j] = 0.0;
        PI1[j] = 1.0;
    }

    nn = 2*nang-1;

    for (int j = 1; j <= nang; j++){
        S1[j] = 0.0;
        S2[j] = 0.0;
    }

    //Riccati-Bessel functins with real argument x calculated by upward recurence
    PSI0=cos(dx);
    PSI1=sin(dx);
    CHI0=-sin(x);
    CHI1=cos(x);
    APSI0=PSI0;
    APSI1=PSI1;
    XI0=APSI0-CHI0*imaginary;
    XI1=APSI1-CHI1*imaginary;
    Qscat_f=0.0;

    int n=1;
    while(n-1-nstop<0){
        DN=n;
        RN=n;
        FN=(2*RN+1)/(RN*(RN+1));
        PSI=(2*DN-1)*PSI1/dx-PSI0;
        APSI=PSI;
        CHI=(2*RN-1)*CHI1/x-CHI0;
        XI=APSI-CHI*imaginary;
        AN=((D[n]/refrel+RN/x)*APSI-APSI1)/((D[n]/refrel+RN/x)*XI-XI1);
        BN=((refrel*D[n]+RN/x)*APSI-APSI1)/((D[n]*refrel+RN/x)*XI-XI1);
        Qscat_f=Qscat_f+(2.0*RN+1.0)*(abs(AN)*abs(AN)+abs(BN)*abs(BN));


        for (int j=1; j<=nang; j++){
            int jj=2*nang-j;
            PI[j]=PI1[j];
            TAU[j]=RN*AMU[j]*PI[j]-(RN+1)*PI0[j];
            P=pow(-1.0,(n-1));
            S1[j]=S1[j]+FN*(AN*PI[j]+BN*TAU[j]);
            T=pow((-1),n);
            S2[j]=S2[j]+FN*(AN*TAU[j]+BN*PI[j]);
            if(j != jj){
                S1[jj]=S1[jj]+FN*(AN*PI[j]*P+BN*TAU[j]*T);
                S2[jj]=S2[jj]+FN*(BN*PI[j]*P+AN*TAU[j]*T);

            }
        }

        PSI0=PSI1;
        PSI1=PSI;
        APSI1=PSI1;
        CHI0=CHI1;
        CHI1=CHI;
        XI1=APSI1-CHI1*imaginary;
        n=n+1;
        RN=n;
        for (int j=1; j<=nang; j++){
            PI1[j]=((2*RN-1)/(RN-1))*AMU[j]*PI[j]-RN*PI0[j]/(RN-1);
            PI0[j]=PI[j];
        }
    }

    Qscat_f=(2.0/(x*x))*Qscat_f;
    Qext_f=(4.0/(x*x))*real(S1[1]);
    Qabs_f = Qext_f-Qscat_f;
    //Note: Qback is not Qbb, but the radar back scattering.
    Qback_f=(4.0/(x*x))*(abs(S1[2*nang-1])*abs(S1[2*nang-1]));
    cout << Qscat_f << endl;
    cout << Qext_f<< endl;
    cout << Qabs_f << endl;

    // Move scattering efficiencies out of the function
    *Qscat_p = Qscat_f; // update the value of Qscat in main using a pointer
    *Qext_p  = Qext_f; // update the value of Qext in main using a pointer
    *Qabs_p = Qabs_f;
    *Qback_p = Qback_f; // update the value of Qback in main using a pointer
    // Move S1 and S2 out of the function
    for (int i=1; i<=nang*2-1; i++){
      S1_p[i] = S1[i];
    }

    for (int i=1; i<=nang*2-1; i++){
      S2_p[i] = S2[i];
    }
    return 0;
}
