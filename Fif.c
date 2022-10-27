/* 
    C implementation of Fast Iterative Filtering (FIF)
    
    written based on  FIF_v2_13.m (MATLAB VERSION)

    Authors: Igor Bertello, Emanuele Papini, Antonio Cicone
    Affiliation(s): IAPS - INAF, University of L'Aquila (Italy)

    Dependencies: FFTW3

*/


#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <fftw3.h>
#include "FFT.h"
#include "Fif.h"
#include "lib/interp.h"


#define DELTA 0.001
#define MaxInner 200
#define Xi 1.6 
#define ExtPoints 3
// the alpha value is set by default to average value (see line 490)
int comp (const void * elem1, const void * elem2) 
{
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}





Maxmins MaxMin;

double *fappo; // vettore di appoggio

//*********************** getMask *******************************
double *getMask(double *y, int n, double k, int *dim_a, double tol)
{
    /*
% Rescale the mask y so that its length becomes 2*k+1.
% k could be an integer or not an integer.
% y is the area under the curve for each bar
*/

    int lim;
    double m;
    double *a, *f_appo;
    double s, s2, t, t1, t2, c;
    double parz, extra, new_k, dx, dy;
    double Norm1;
    bool f_appoB = 0;

    //n=length(y);
    m = (double)(n - 1) / 2.;

    if (k <= m)
    { //% The prefixed filter contains enough points

        if ((int)k == k)
        { //  % if the mask_length is an integer
            lim = 2 * (int)k + 1;

            a = (double *)malloc(sizeof(double) * (lim)); //zeros(1,2*k+1);

            for (int i = 1; i < lim+1; i++)
            {
                s = ((double)(i - 1) * (2 * m + 1)) / (double)(lim) + 1;
                t = (double)(i) * (2 * m + 1)/ (double)(lim);
                s2 = ceil(s) - s;
                t1 = t - floor(t);

                if (floor(t) < 1)
                {
                    printf("#Ops\n");
                }
                parz = 0;
                for (int j = (int)ceil(s) - 1; j < (int)floor(t); j++)
                { //    a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
                    parz += y[j];
                }
                a[i-1] = parz + s2 * y[(int)ceil(s)-1] + t1 * y[(int)floor(t)-1];
            }
            *dim_a = lim;
        }
        else
        { // % if the mask length is not an integer
            new_k = floor(k);
            extra = k - new_k;
            c = (2 * m + 1) / (2 * new_k + 1 + 2 * extra);

            a = (double *)malloc(sizeof(double) * (2 * (int)new_k + 3)); //zeros(1,2*new_k+3);

            t = extra * c + 1;
            t1 = t - floor(t);
            //%t2=ceil(t)-t;
            if (k < 0)
            {
                // disp('Ops')
                free(a); //a=[];
                return 0;
            } //end
            parz = 0;
            for (int i = 0; i < (int)floor(t); i++)
            { //    a(1)=sum(y(1:floor(t)))+t1*y(floor(t));
                parz += y[i];
            }
            a[0] = parz + t1 * y[(int)floor(t)-1];

            for (int i = 2; i < (int)(2 * new_k) + 3; i++)
            {
                s = extra * c + (double)(i - 2) * c + 1;
                t = extra * c + (double)(i - 1) * c;
                s2 = ceil(s) - s;
                t1 = t - floor(t);

                parz = 0;
                for (int i = (int)ceil(s)-1; i < (int)floor(t); i++)
                { //  a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
                    parz += y[i];
                }
                a[i-1] = parz + s2 * y[(int)ceil(s)-1] + t1 * y[(int)floor(t)-1];
            }
            t2 = ceil(t) - t;

            parz = 0;
            for (int i = (int)ceil(t)-1; i < n; i++)
            { //   a(2*new_k+3)=sum(y(ceil(t):n))+t2*y(ceil(t));
                parz += y[i];
            }
            a[(int)(2 * new_k) + 2] = parz + t2 * y[(int)ceil(t)-1];
            *dim_a = 2 * (int)new_k + 3;//(int)floor(t);
        } // mod(k,1) end
    }
    else
    { 

        dx = 0.01;
        f_appo = (double *)malloc(sizeof(double) * n);
        f_appoB = 1;
        
        double *xappo =(double *)malloc(sizeof(double) * ((int)m+1));
        for (int i = 0; i < ((int)m+1); i++){xappo[i] = (double)i;}

        for (int i = 0; i < n; i++){f_appo[i] = y[i] / dx;}
        dy = m * dx / k;
        
        
        lim = (int)floor(k) ; 
        double *xappo_interp =(double *)malloc(sizeof(double) * (lim));
        for (int i = 0; i < lim; i++){xappo_interp[i] =(double)i*m/k;}
        
        double *f_appo_interp;         
        f_appo_interp = interp_linear ( 1, (int)m+1, xappo, &f_appo[(int)m],lim, xappo_interp);

        a = (double *)malloc(sizeof(double) * (2 * lim -1)); //zeros(1,2*new_k+3);
        *dim_a = 2*lim -1; 

        for (int i = 0; i < lim-1; i++){a[i] = f_appo_interp[lim-1-i];}
        memcpy(a+lim-1, f_appo_interp, lim * 8); // copia la maschera traslata di 
        free(f_appo_interp); free(xappo_interp); free(xappo);

        Norm1 = 0; //norm(f,1);
        for (int i = 0; i < *dim_a; i++)
        {
            Norm1 += abs(a[i]);
        }

        if (abs((Norm1) - 1) > tol)
        {
            //if verbose>0  {
            printf("\n\n Warning!\n\n");
            printf(" Area under the mask equals %2.20f\n", Norm1);
            printf(" it should be equal to 1\n We rescale it using its norm 1\n\n");
            // } //end
            for (int i = 0; i < *dim_a; i++)
                a[i] = a[i] / Norm1;
        } //end*/
    }     // end*/

    if (f_appoB)
        free(f_appo);
    return a;
}

//********************** Maxmins_v3_8 *******************************

Maxmins Maxmins_v3_8(double *f, unsigned int N, double tol)
{

    /* v3_8 Based on version Maxmins_v3_8.m
 Identify the maxima and minima of a signal f


*/

    unsigned int *Maxs; // i massimi, dimensionamento dinamico
    unsigned int *Mins; // i minimi, dimensionamento dinamico
    double *df;         // il vettore differenze, dimensionamento dinamico
    double *f1, *z;     // Funzione di appoggio, dimensione dinamica
    unsigned int N_old, last_df;
    double fappo1, fappo2, fappo3; // float di appoggio
    int h, cmaxs, cmins, c;
    int posc;                   // posizione del max/min
    int ctot, mac = 0, mic = 0; // indici per i vettori Maxs, Mins
    Maxmins valori;
    bool dfB = 0, f1B = 0;
    bool MMB = 0, MmB = 0;

    // dimensionamento vettori
    Maxs = (unsigned int *)malloc(sizeof(unsigned int) * N);
    Mins = (unsigned int *)malloc(sizeof(unsigned int) * N);

    df = (double *)malloc(sizeof(double) * (N-1));

    MMB = MmB = dfB = 1;

    //Creo il vettore differenza

    for (int x = 0; x < N - 1; x++)
    {
        df[x] = f[x + 1] - f[x];
    }

    h = 1;

    while ((h < N) && (abs(df[h-1] / f[h-1]) <= tol))
        h++; // trova il primo punto non costante.  df[x]/f[x]

    if (dfB)
    {
        free(df);
        dfB = 0;
    }


    cmaxs = 0;
    cmins = 0;
    c = 0;

    N_old = N;

    if ((f1 = (double *)malloc(sizeof(double) * (N + h ))) != NULL)
        f1B = 1; //f=[f f(2:h)];
    if ((df = (double *)malloc(sizeof(double) * (N + h ))) != NULL)
        dfB = 1;                              //df=diff([f f(2:h+1)]);

    for (int i = 0; i < N; i++)
    { //memcpy(f1, f, N*8);
        f1[i] = f[i];
    }
    for (int i = 0; i < h; i++)
    { //memcpy(f1+(N*8), fappo, h*8);  // copio in fondo a F1 la parte costante di lughezza h
        f1[N + i] = f[i+1];
    }
    N = N + h; // +1?

    //ricreo il vettore differenza
    for (int x = 0; x < N - 1; x++)
    {
        df[x] = f1[x + 1] - f1[x];
    }

    for (int i = h-1; i < N - 2; i++)
    { 

        fappo2 = ((df[i] * df[i + 1]) / (abs(f1[i]) * abs(f1[i]))); // differenza relativa al quadrato: ~derivata seconda

        if ((fappo2 <= tol) && (fappo2 >= -tol))
        {                                //  ~Derivata seconda=0
            fappo3 = df[i] / abs(f1[i]); // diffferenza relativa i-sima: ~derivata prima
            if (fappo3 < -tol)
            {
                last_df = -1;
                posc = i;
            }
            else if (fappo3 > tol)
            {
                last_df = +1;
                posc = i;
            }
            else if (df[i] == 0)
            {
                last_df = 0;
                posc = i;
            }
            c++;
            fappo1 = df[i + 1] / abs(f1[i]); // diffferenza relativa (i+1)-sima: ~derivata prima
            if (fappo1 < -tol)
            {
                if ((last_df == +1) || (last_df == 0))
                {
                    Maxs[cmaxs] = ((posc + (int)(floor((double)(c - 1) / 2)) + 1) % (N_old-1));
                    cmaxs++; //cmaxs=cmaxs+1;
                }
                c = 0;
            }
            if (fappo1 > tol)
            {
                if ((last_df == -1) || (last_df == 0))
                {
                    Mins[cmins] = ((posc + (int)(floor((double)(c - 1) / 2)) + 1) % (N_old-1));
                    cmins++; //cmins=cmins+1;
                }
                c = 0;
            }

        } // entro tol: ~0
        if (fappo2 < -tol)
        {
            fappo3 = df[i] / abs(f1[i]);     // diffferenza relativa i-sima: ~derivata prima
            fappo1 = df[i + 1] / abs(f1[i]); // diffferenza relativa (i+1)-sima: ~derivata prima

            if ((fappo3 < -tol) && (fappo1 > tol))
            {

                Mins[cmins] = ((i + 1) % (N_old-1));
                cmins++;
                last_df = -1;
            }
            else if ((fappo3 > tol) && (fappo1 < -tol))
            {
                Maxs[cmaxs] = ((i + 1) % (N_old-1));
                cmaxs++;
                last_df = +1;
            }
        }
    } // end for i*/

    if (c > 0)
    {

        if ((cmins > 0) && (Mins[cmins] == 0))
        {
            Mins[cmins] = N;
        }
        if ((cmaxs > 0) && (Maxs[cmaxs] == 0))
        {
            Maxs[cmaxs] = N;
        }
    } // end

    // Crea il vettore di output maxmins
    ctot = cmaxs + cmins;

    valori.nout = ctot;
    valori.maxmins = (unsigned int *)malloc(sizeof(unsigned int) * ctot);
    
    for (int i=0;i<cmaxs;i++){valori.maxmins[i] = Maxs[i];}
    for (int i=0;i<cmins;i++){valori.maxmins[i+cmaxs] = Mins[i];}
    qsort(valori.maxmins, ctot, sizeof(unsigned int), comp);

    // libera la memoria
    if (f1B)
    {
        free(f1);
        f1B = 0;
    }
    if (dfB)
    {
        free(df);
        dfB = 0;
    }
    if (MmB)
    {
        free(Mins);
        MmB = 0;
    }
    if (MMB)
    {
        free(Maxs);
        MMB = 0;
    }
    free(z);

    return valori;
} // end function

//********************** FIF_v2_1 *******************************
Fif_t FIF_v2_1(double *f, int N, int *maxIMF)
{

    int numIMF, dim_a, nh;
    unsigned int k_pp;
    int in_step, ext_sig;
    int posF; 
    double df;
    int Nza, Nxs, N_old, N_r;
    int cont;
    double rho, theta, thetan;
    double Norm1 = 1, m, logM = 0;
    double SD, norm_n, norm_d; // numeratore e denominatore per il calolo di SD
    double *h;
    double *mask, *mappo;
    double *ifftA, *fappo, *z;
    int n ;

    n = *maxIMF;


    FILE *fpin;
    unsigned int conta = 0, size;
    bool hB, ifftAB, fappoB, maskB, mappoB;

    fftwl_complex *fftH, *ffth_new, *ffth_old;
    Fif_t IMF, *IMFappo, *next;
    double tol = 10e-12;
    hB = ifftAB = fappoB = maskB = mappoB = 0;


    IMFappo = &IMF;
    size = sizeM;

    // / normalizzazione e verifica grandezza del segnale
    Norm1 = 0; // MATLAB:  norm(f,1);
    for (int i = 0; i < N; i++)
    {
        Norm1 += abs(f[i]);
    }
    for (int i = 0; i < N; i++)
    {
        f[i] /= Norm1;
    }

    N_r = 0;
    fappo = (double *)malloc(sizeof(double) * (N + 11));
    fappoB = 1;
    for (int i = 0; i < N; i++)
    { //   memcpy(fappo, f, n*8);
        if (abs(f[i]) > tol)
        { // copia f solo se il valore è magigore di tol
            fappo[N_r++] = f[i];
        }
    }
    if (N_r == 0)
        return IMF; // se l'intero segnale è troppo debole esce
    for (int i = 0; i < 10; i++)
        fappo[N_r + i] = fappo[i]; //memcpy(fappo+N_r,fappo, 80);


    MaxMin = Maxmins_v3_8(fappo, N_r + 10, tol);

    if (fappoB)
    {
        free(fappo);
        fappoB = 0;
    }
    k_pp = MaxMin.nout;
    if (k_pp == 0)
    {

        return IMF;
    }
    for (int i=0;i<MaxMin.nout;i++){
    printf("%d ",MaxMin.maxmins[i]);
    }
    printf("\n");

    numIMF = 0;
    cont = 0; // numero di IMF fatte
    
    while ((numIMF < n) && (k_pp > ExtPoints))
    { // 1. ciclo su tutte le IMF:   numIMF<n

        numIMF++; // siamo alla IMF successiva
        cont++;
        SD = 1;
        h = (double *)malloc(sizeof(double) * N);
        hB = 1;
        memcpy(h, f, sizeof(double) * N);
        nh = N;
        
        //for (int i = 0; i < N; i++){printf("H %d, %f\n", i, h[i]);}
        // 2. definizione lunghezza maschera

        m =  round(2*N_r / (k_pp + 1) * Xi); // MATLAB:  2*round(N_pp/k_pp*options.Xi);

        if (numIMF > 1)
        {   // 2.a  c'è almeno una IMF fatta
            // / MATLAB:  logM=stats(countIMFs-1).logM
            if (m <= logM)
            { // MATLAB:  m<=stats(countIMFs-1).logM).

                if (1)
                { // 2.c options.MonotoneMaskLength==true
                    m = ceil(logM * 1.1);
                } // 2.c  opzione monotoneMaskLength
            }     // 2.b MATLAB:  m<logM(i-1)
        }         // 2.a almeno una IMF fatta

        logM = m;
        IMFappo->stats.logM = m;
        printf("%f \n",IMFappo->stats.logM);//ema
        /// ****************
        // crea la maschera ridotta mask
        if (maskB)
        {
            maskB = 0;
            free(mask);
        }
        mask = (double *)malloc(sizeof(double) * N);
        mask = getMask(MM, size, m, &dim_a, tol); // MATLAB:   a = get_mask_v1_1(MM,m,options.verbose,tol);
        mask = (double *)realloc(mask, sizeof(double) * dim_a);

        ext_sig = 0;
        N_r = N;
        
        if (N < dim_a)
        { //if (N<length(a)){  // 3. estensione segnale
            ext_sig = 1;
            Nxs = (int)ceil((float)dim_a / (float)N);
            N_old = N;
            if ((Nxs % 2) == 0)
            { // Nxs è pari
                Nxs++;
            }
            N_r = N * Nxs;

            if (fappoB)
            {
                free(fappo);
                fappoB = 0;
            }
            fappo = (double *)malloc(sizeof(double) * N_r);
            fappoB = 1;

            for (int i = 0; i < Nxs; i++)
            {
                memcpy(fappo + N * i, h, N * 8);
            }
            if (hB)
            {
                free(h);
                hB = 0;
            }
            h = (double *)malloc(sizeof(double) * N_r);
            hB = 1;
            memcpy(h, fappo, 8 * N_r);

            if (fappoB)
            {
                free(fappo);
                fappoB = 0;
            }
        } // 3. fine if estensione segnale
        Nza = N_r - dim_a; // attenzione se viene negativo!!
        printf("N_r: %d, dim_a: %d, Nza : %d\n",N_r,dim_a,Nza);
        // / Fa la FFT della maschera mask
        if (mappoB)
        {
            free(mappo);
            mappoB = 0;
        }
        mappo = (double *)malloc(sizeof(double) * N_r); // a = [zeros(1,(Nza-1)/2) a zeros(1,(Nza-1)/2+1)];
        mappoB = 1;

        for (int i = 0; i < N_r; i++){mappo[i] = 0;} // mette a zero 
        if ((Nza % 2) == 0)
        { // 4a. Nza è pari. iffT
            //******* MATLAB:  ifftA=real(fft([a((length(a)-1)/2+1:end) a(1:(length(a)-1)/2)]));
            memcpy(mappo, mask + dim_a/2, dim_a/2 * 8); // copia la maschera traslata di 
            memcpy(mappo+(N_r - dim_a/2 ),mask,dim_a/2 * 8);
        } // 4a. fine Nza è pari. iffT
        else
        { // 4b. Nza è dispari. iffT
            memcpy(mappo, mask + dim_a/2, (dim_a/2+1) * 8); // copia la maschera traslata di 
            memcpy(mappo+(N_r - dim_a/2 ),mask,dim_a/2 * 8);

        } // 4b. fine Nza è dispari. iffT

        if (ifftAB)
        {
            free(ifftA);
            ifftAB = 0;
        }
        ifftAB = 1;
        fflush(NULL);
        ifftA = realFFT(mappo, N_r); // max 121. 122 crash
        fflush(NULL);
        if (mappoB)
        {
            free(mappo);
            mappoB = 0;
        }

        // / fa la FFT del segnale
        ffth_new = (fftwl_complex *)fftwl_malloc(sizeof(fftwl_complex) * N_r);
        ffth_old = (fftwl_complex *)fftwl_malloc(sizeof(fftwl_complex) * N_r);
        fftH = fft_dir(h, N_r); // faccio la FFT del segnale h

        //MATLAB: posF=find((diff(ifftA)>0)~=0,1,'first');
        posF = 0;
        df = -1.;
        while ((posF < N_r - 1) && (df < 0))
        {
            df = ifftA[posF + 1] - ifftA[posF];
            posF++;
        }
        posF--;
        printf("Posf : %d\n",posF);

        IMFappo->stats.posF = posF;
        IMFappo->stats.valF = ifftA[posF];

        for (int i = 0; i < N_r; i++)
        { //   ifftA=ifftA-ifftA(posF);
            ifftA[i] -= IMFappo->stats.valF;
            if (ifftA[i] < 0)
                ifftA[i] = 0; //  ifftA(ifftA<0)=0;
        }
        
        in_step = 0;
        while ((SD > DELTA) && (in_step < MaxInner))
        {                 // 5.   nucleo dell'algoritmo. delta=0.001 MaxInner=200
            in_step += 1; // diverso da +1 serve a far convergere il calcolo più in fretta

            // / fa la convoluzione tra ifftA e fftH

            norm_n = norm_d = 0;
            for (int i = 0; i < N_r; i++)
            {
                ffth_new[i][0] = pow((1 - ifftA[i]), (in_step)) * fftH[i][0];
                ffth_new[i][1] = pow((1 - ifftA[i]), (in_step)) * fftH[i][1];

                ffth_old[i][0] = pow((1 - ifftA[i]), (in_step - 1)) * fftH[i][0];
                ffth_old[i][1] = pow((1 - ifftA[i]), (in_step - 1)) * fftH[i][1];
                norm_n += (ffth_new[i][0] - ffth_old[i][0]) * (ffth_new[i][0] - ffth_old[i][0]) + (ffth_new[i][1] - ffth_old[i][1]) * (ffth_new[i][1] - ffth_old[i][1]);
                norm_d += (ffth_old[i][0] * ffth_old[i][0]) + (ffth_old[i][1] * ffth_old[i][1]);
            }

            //       %%%%%%%%%%%%%%%% Updating stopping criterium %%%%%%%%%%%%%%%%%
            SD = norm_n / norm_d; // MATLAB: SD=norm(fft_h_new-fft_h_old)^2/norm(fft_h_old)^2;
            //printf("step #: %d ,SD: %lf \n",in_step,SD);
        } // 5. fine nucleo dell'algoritmo

        // / fa la trasformata inversa di fft_h_new
        if (hB)
        {
            free(h);
            hB = 0;
        }
        //h = (double *)malloc(sizeof(double) * N_r);

        h = fft_inv(ffth_new, N_r);
        hB = 1;
        if (ext_sig)
        { // % we reduce the signal
            N = N_old;
            // fappo=(double *)malloc(sizeof(double)*N);
            int j = 0;
            for (int i = (N * (Nxs - 1) / 2 ); i < (N * ((Nxs - 1) / 2 + 1)); i++)
            { // /  h=h(N*(Nxs-1)/2+1:N*((Nxs-1)/2+1));
                h[j++] = h[i];
            }
            if (j == N)
                h = (double *)realloc(h, sizeof(double) * N);
            else
                printf("Error in shrinking the signal from the extension!\n");
                
        }

        IMFappo->stats.in_step = in_step;


        IMFappo->dati = (double *)malloc(sizeof(double) * N); 
        //copying result and renormalizing to the original amplitude
        for (int i=0;i<N;i++){IMFappo->dati[i] = h[i]*Norm1;}
        
        for (int i = 0; i < N; i++)
        { // f=f-h;
            f[i] -= h[i];
        }
        next = (Fif_t *)malloc(sizeof(Fif_t)); // aggiunge un elemnto alla lista delle IMF
        IMFappo->next = next;           
        IMFappo = next;
        N_r = 0;

        fappo = (double *)malloc(sizeof(double) * N + 80);
        fappoB = 1;

        for (int i = 0; i < N; i++)
        { //   memcpy(fappo, f, n*8);
            if (abs(f[i]) > tol)
            { // copia f solo se il valore è maggiore di tol
                fappo[N_r++] = f[i];
            }
        }

        if (N_r == 0)
            return IMF; // se l'intero segnale è troppo debole esce
        memcpy(fappo + N_r, fappo, 80);

        if (logM >= 20)
        {
            // deve fare la media mobile
            double *sappo = (double *)malloc(sizeof(double) * N + 80);
            sappo[0] = fappo[0]; 
            sappo[1] = (fappo[0] + fappo[1])/2;
            sappo[2] = (fappo[0]+ fappo[1] +fappo[2]+fappo[3])/4;
            sappo[3] = (fappo[0]+ fappo[1] +fappo[2]+fappo[3] +fappo[4]+fappo[5])/6;
            sappo[4] = (fappo[0]+ fappo[1] +fappo[2]+fappo[3]+fappo[4] + 
                        fappo[5] +fappo[6]+fappo[7])/8;
            sappo[N_r+ 9] = fappo[N_r+9]; 
            sappo[N_r+ 8] = (fappo[N_r+9] + fappo[N_r+8])/2;
            sappo[N_r+ 7] = (fappo[N_r+9] + fappo[N_r+8] +fappo[N_r+7]+fappo[N_r+6])/4;
            sappo[N_r+ 6] = (fappo[N_r+9] + fappo[N_r+8] +fappo[N_r+7]+fappo[N_r+6] + 
                             fappo[N_r+5]+fappo[N_r+4])/6;
            sappo[N_r+ 5] = (fappo[N_r+9] + fappo[N_r+8] +fappo[N_r+7]+fappo[N_r+6] + 
                             fappo[N_r+5]+fappo[N_r+4]+fappo[N_r+3]+fappo[N_r+2])/8;
            for (int i = 5; i < N_r + 5; i++)
            {
                sappo[i] = (fappo[i-5] + fappo[i-4] + fappo[i-3] + fappo[i-2] + fappo[i-1] +
                           fappo[i] + fappo[i+1] + fappo[i+2] + fappo[i+3] + fappo[i+4])/10;
            }
            memcpy(fappo,sappo, 8*N_r+80);
            MaxMin = Maxmins_v3_8(fappo, N_r + 10, tol); //  maxmins_pp=Maxmins_v3_8(movmean([f_pp f_pp(1:10)],10),tol);
            free(sappo);
        }
        else
        {
            MaxMin = Maxmins_v3_8(fappo, N_r + 80, tol);
        }

        k_pp = MaxMin.nout;

        if (k_pp == 0)
        {
            printf("#No extrema detected. numIMF=%d k_pp=%d N_r=%d\n", numIMF, k_pp, N_r);
            return IMF;
        }
    } 
    //saving residual of the decomposition in the last imf

    next = (Fif_t *)malloc(sizeof(Fif_t)); // aggiunge un elemnto alla lista delle IMF
    IMFappo->next = next;           
    IMFappo = next;
    for (int i=0;i<N;i++){IMFappo->dati[i] = f[i]*Norm1;}

    *maxIMF = numIMF;
    fflush(NULL);
    return IMF;
}
