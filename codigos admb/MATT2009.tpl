GLOBALS_SECTION
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc.csv");

TOP_OF_MAIN_SECTION
 time(&start);
 arrmblsize = 90000000;
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
 gradient_structure::set_MAX_NVAR_OFFSET(5000);
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


DATA_SECTION
 init_int ntime
 init_int nedades
 init_int ntallas
 init_int nper

 init_matrix mdatos(1,ntime,1,13) // Yrs Captura cv1 CPUEI cv2 CPUEA cv3 Bcru cv4 MPH cv5 nmf nmc
 init_vector Tallas(1,ntallas)
 init_matrix Ctot(1,ntime,1,ntallas)
 init_matrix Ncru(1,ntime,1,ntallas)
 init_vector msex(1,ntallas)
 init_vector Wmed(1,ntallas)

//!!cout<<Wmed<<endl;exit(1);
//!! ad_comm::change_datafile_name("ModtallasLL.ctl");
 init_number sigmaR
 init_vector dt(1,3)
 init_vector Par_bio(1,7)    // Loo k Lr sr b M h
 init_vector cv_Par_bio(1,7) // soo activo Lr
 init_number priorlog_Rmed
  //!!cout<<cv_Par_bio<<endl;exit(1);

 init_int    minedad //edad madurez
 init_number bprior  //coef Hiperestabilidad/reduccion (1=proporcional)

  number log_Lr_prior
  number log_sr_prior
  number log_b_prior
  number log_beta_prior

  !! log_Lr_prior   = log(Par_bio(3));
  !! log_sr_prior   = log(Par_bio(4));
  !! log_beta_prior = log(Par_bio(5));
  !! log_b_prior    = log(bprior);

 init_number L50prior  // (L50= talla modal	app)
 init_number s1prior
 init_number s2prior  // (modelo domo)
 init_int opt_sel2    // opcion domo en flota
 init_int opt_sel3    // opcion domo en crucero (<0, logistico, >0, domo)

 number log_L50prior
 number log_s1prior
 number log_s2prior

 !! log_L50prior = log(L50prior);
 !! log_s1prior = log(s1prior);
 !! log_s2prior = log(s2prior);

 init_int    nbloques1
 init_vector ybloques1(1,nbloques1)

 init_int    nbloques2
 init_vector ybloques2(1,nbloques2)

 init_int    nqbloques
 init_vector yqbloques(1,nqbloques)

 init_int    nqbloquesc
 init_vector yqbloquesc(1,nqbloquesc)

 init_int    nqbloquesmph
 init_vector yqbloquesmph(1,nqbloquesmph)


 init_number      prior_qfI
 init_number      prior_qfA
 init_number      prior_qc
 init_number      prior_qmph

 number      log_prior_qfI
 number      log_prior_qfA
 number      log_prior_qc
 number      log_prior_qmph

 !!log_prior_qfI=log(prior_qfI);
 !!log_prior_qfA=log(prior_qfA);
 !!log_prior_qc=log(prior_qc);
 !!log_prior_qmph =log(prior_qmph);

 init_vector cv_q(1,4)

 init_int    opt_qfI
 init_int    opt_qfA
 init_int    opt_qc
 init_int    opt_qmph

 init_int    opt1_fase // selectividad flota
 init_int    opt2_fase // selectividad cruceros


 init_int    opt_Lr
 init_int    opt_sr
 init_int    opt_beta
 init_int    opt_F
 init_int    opt_Ro
 init_int    opt_devRt
 init_int    opt_devNo//Condicion inicial (Si no estima (<0) = poblaci�nen equilibrio)
 init_int    opt_bpow //hiperestabilidad deber�an ser por flota?

 init_int    npbr
 init_vector pbr(1,npbr)
 init_int    ntime_sim

 init_number festim //Opcion_estima_(1)_o_calcula_(0)


INITIALIZATION_SECTION
  log_Rmed       9.43
  log_Rmed       8.4
  log_Lr         log_Lr_prior
  log_sr         log_sr_prior
  log_L50        log_L50prior
  log_L50c       log_L50prior
  log_sigma1     log_s1prior
  log_sigma2     log_s2prior
  log_sigma1c    log_s1prior
  log_sigma2c    log_s2prior
  log_b          log_b_prior
  log_beta       log_beta_prior

  log_qfloI      log_prior_qfI
  log_qfloA      log_prior_qfA
  log_qcru       log_prior_qc
  log_qmph       log_prior_qmph


PARAMETER_SECTION

// selectividad param�trica a la talla com�n
// init_bounded_vector log_L50f(1,nbloques1,-5,8,opt1_fase)

 init_vector log_L50(1,nbloques1,opt1_fase)
 init_vector log_sigma1(1,nbloques1,opt1_fase)
 init_vector log_sigma2(1,nbloques1,opt_sel2)


 init_vector log_L50c(1,nbloques2,opt2_fase)
 init_vector log_sigma1c(1,nbloques2,opt2_fase)
 init_vector log_sigma2c(1,nbloques2,opt_sel3)

// parametros reclutamientos y mortalidades)
 init_number log_Rmed(opt_Ro)
 init_number log_Rmed_pre(1)
 init_bounded_dev_vector log_desv_Rt(1,ntime,-10,10,opt_devRt)
 init_bounded_vector log_desv_No(1,nper,-10,10,opt_devNo)
 init_bounded_vector log_F(1,ntime,-20,0.7,opt_F) // log  mortalidad por pesca por flota

// capturabilidades
 init_vector log_qfloI(1,nqbloques,opt_qfI)
 init_vector log_qfloA(1,nqbloques,opt_qfA)
 init_vector log_qcru(1,nqbloquesc,opt_qc)
 init_vector log_qmph(1,nqbloquesmph,opt_qmph)

 init_number log_b(opt_bpow)

// Crecimiento
 init_number log_Lr(opt_Lr)
 init_number log_sr(opt_sr)
 init_number log_beta(opt_beta)

//---------------------------------------------------------------------------------
//Defino las variables de estado
 vector BMflo(1,ntime)
 vector BMcru(1,ntime)
 vector Brec(1,ntime)
 vector BMmph(1,ntime)
 sdreport_vector pred_CPUEI(1,ntime);
 sdreport_vector pred_CPUEA(1,ntime);
 sdreport_vector pred_Bcru(1,ntime);
 sdreport_vector pred_Desemb(1,ntime);
 vector likeval(1,15);
 vector Neq(1,ntallas);

 sdreport_vector Rpred(1,ntime);
 sdreport_number Rnpre;
 sdreport_number Rpop;

 vector Unos_edad(1,nedades);
 vector Unos_year(1,ntime);
 vector Unos_tallas(1,ntallas);
 vector delta(1,ntallas)
 vector Lesp(1,ntallas)
 vector sigmaL(1,ntallas)
 vector pre(1,ntallas)

 vector mu_edad(1,nedades)
 vector sigma_edad(1,nedades)
 vector BDo(1,ntime);
 vector No(1,ntallas)
 vector prior(1,7)
 vector yrs(1,ntime)
 vector Desemb(1,ntime);
 vector CPUEI(1,ntime);
 vector CPUEA(1,ntime);
 vector Bcru(1,ntime);
 vector mph(1,ntime);
 vector Lmed_obs(1,ntime)
 sdreport_vector Lmed_pred(1,ntime)
 vector Lmed_obsc(1,ntime)
 sdreport_vector Lmed_predc(1,ntime)
 vector edades(1,nedades)
 sdreport_vector Reclutas(1,ntime)
 sdreport_vector Reclut_R(1,ntime)
 vector nm(1,ntime)
 vector nmc(1,ntime)
 vector penalty(1,6)

 matrix cv_index(1,5,1,ntime)

 matrix S1(1,nbloques1,1,ntallas)
 matrix S2(1,nbloques2,1,ntallas)

 matrix Sel(1,ntime,1,ntallas)
 matrix Selc(1,ntime,1,ntallas)

 matrix F(1,ntime,1,ntallas)
 matrix Z(1,ntime,1,ntallas)
 matrix S(1,ntime,1,ntallas)


 matrix N(1,ntime,1,ntallas)

 matrix NM(1,ntime,1,ntallas)
 matrix NMD(1,ntime,1,ntallas)
 matrix NDv(1,ntime,1,ntallas)
 matrix Nrec(1,ntime,1,ntallas)
 matrix NVflo(1,ntime,1,ntallas)
 matrix NVcru(1,ntime,1,ntallas)
 matrix NVmph(1,ntime,1,ntallas)

 matrix pred_Ctot(1,ntime,1,ntallas)

 matrix pobs(1,ntime,1,ntallas)
 matrix ppred(1,ntime,1,ntallas)
 matrix pobsc(1,ntime,1,ntallas)
 matrix ppredc(1,ntime,1,ntallas)

 matrix T(1,ntallas,1,ntallas)

 matrix Nv(1,ntime,1,nedades)
 matrix NMDv(1,ntime,1,ntallas)

 number suma1
 number suma2
 number suma3
 number suma4
 number suma5
 number suma6

 number So
 number alfa
 number beta

 number Linf
 number k
 number Linfh
 number M
 number Lr
 number sr
 number Lm
 number Rm
 number h

 number BDp
 number Npplus
 number Bp_anch

 number nm1;
 number cuenta1;
 number alfa_sr;
 number beta_sr;
 number pF

 vector Np(1,ntallas)
 vector Zpbr(1,ntallas)
 vector Fpbr(1,ntallas)
 vector Sp(1,ntallas)

 matrix Bp(1,npbr,1,ntime_sim)
 vector CTPp(1,ntallas)
 matrix Rpp(1,npbr,1,ntime_sim)


 objective_function_value f

 sdreport_vector BD(1,ntime) //
 sdreport_vector BT(1,ntime) //
 sdreport_vector RPRlp(1,ntime) //
 sdreport_vector pred_mph(1,ntime);
 sdreport_number SSBo
 sdreport_matrix Yp(1,npbr,1,ntime_sim)
 sdreport_vector Ypact(1,npbr)
 sdreport_vector RPRequ3(1,ntime)
 sdreport_vector Frpr(1,ntime)

PRELIMINARY_CALCS_SECTION

 yrs=column(mdatos,1);
 Desemb=column(mdatos,2);
 CPUEI=column(mdatos,4);
 CPUEA=column(mdatos,6);
 Bcru=column(mdatos,8);
 mph=column(mdatos,10);
 nm=column(mdatos,12);
 nmc=column (mdatos,13);

 edades.fill_seqadd(minedad,1);

 cv_index(1)=column(mdatos,3);
 cv_index(2)=column(mdatos,5);
 cv_index(3)=column(mdatos,7);
 cv_index(4)=column(mdatos,9);
 cv_index(5)=column(mdatos,11);

 Linf=Par_bio(1);
 k=Par_bio(2);
 M=Par_bio(6);
 h=Par_bio(7);

 Unos_tallas=1; // lo uso en operaciones matriciales con tallas
 Unos_year=1;   // lo uso en operaciones matriciales con el a�o


RUNTIME_SECTION
  convergence_criteria 1.e-1,1.e-01,1.e-03,1e-3,1e-5
  maximum_function_evaluations 100,100,200,3000,3500

PROCEDURE_SECTION
// se listan las funciones que contienen los calculos
 Eval_Trans_talla_talla();
 Eval_selectividad();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 Eval_logverosim();
 Eval_funcion_objetivo();

 if(last_phase()){Eval_CTP();}


//-----------------------------------------------------------------
FUNCTION Eval_Trans_talla_talla

  Linf=Par_bio(1);
  k=Par_bio(2);
  beta=Par_bio(5);

//  if(active(log_k)){k=mfexp(log_k);}
  if(active(log_beta)){beta=mfexp(log_beta);}

 int i, j;

// matriz de transicion modelo normal

  delta=(Linf-Tallas)*(1-mfexp(-k));// incremento en tallas
  Lesp=Tallas+delta; // talla esperada luego del crecimiento
  sigmaL=delta*beta;

  for (i=1;i<=ntallas;i++){
    for (j=1;j<=ntallas;j++){
      if(i==j){
         T(i,j)=1.0;}}
   }


  for (i=1;i<=ntallas;i++){

    for (j=1;j<=ntallas;j++){
     if(sigmaL(i)>0){
     T(i,j)=mfexp(-0.5*square((Lesp(i)-Tallas(j))/sigmaL(i)));}}
   }


  for (j=1;j<=ntallas;j++){
  T(j)/=sum(T(j));
  }


//----------------------------------------------------------------------

FUNCTION Eval_selectividad
 int i,j;

 // FLOTA...................

 for (j=1;j<=nbloques1;j++){

 S1(j)=exp(-0.5*square(Tallas-exp(log_L50(j)))/square(exp(log_sigma1(j))));


    for (i=1;i<=ntallas;i++){

      if(Tallas(i)>=exp(log_L50(j))){
      S1(j,i)= exp(-0.5*square(Tallas(i)-exp(log_L50(j)))/square(exp(log_sigma2(j))));
      }

 }}

   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques1;j++){
              if (yrs(i)>=ybloques1(j)){
                Sel(i)=S1(j);}
       }
   }

 // CRUCERO...................

 for (j=1;j<=nbloques2;j++){

 S2(j)=exp(-0.5*square(Tallas-exp(log_L50c(j)))/square(exp(log_sigma1c(j))));

    for (i=1;i<=ntallas;i++){

      if(Tallas(i)>=exp(log_L50c(j))){
      S2(j,i)= exp(-0.5*square(Tallas(i)-exp(log_L50c(j)))/square(exp(log_sigma2c(j))));
      }

 }}


   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques2;j++){
              if (yrs(i)>=ybloques2(j)){
                Selc(i)=S2(j);}
       }
   }


FUNCTION Eval_mortalidades

 F=elem_prod(Sel,outer_prod(mfexp(log_F),Unos_tallas));

 Z=F+M;

 S=mfexp(-1.0*Z);


FUNCTION Eval_abundancia
 int i, j;

  Lr=Par_bio(3);
  sr=Par_bio(4);

  if (active(log_Lr)){Lr=mfexp(log_Lr);}
  if (active(log_sr)){sr=mfexp(log_sr);}


// genero la composicion de tallas del reclutamiento
  pre=exp(-0.5*square((Tallas-Lr)/sr));
  pre/=sum(pre);

// Para no estimar Rmed
  if(opt_Ro<0){
  log_Rmed=priorlog_Rmed;
  }

// genero una estructura inicial en torno a Z del primer a�o;
  Reclut_R = mfexp(log_Rmed+log_desv_Rt);
  Reclutas=mfexp(log_Rmed+log_desv_Rt);

// genero la poblacion en equilibrio virginal de LP;

  No=pre*exp(log_Rmed);
  for (int j=1;j<=nper*nedades;j++){
  No=(No*exp(-1.*M))*T+pre*exp(log_Rmed);
  }
  SSBo=sum(elem_prod(No*mfexp(-dt(1)*M),elem_prod(Wmed,msex)));
  alfa_sr=4*h*exp(log_Rmed+0.5*square(sigmaR))/(5*h-1);//
  beta_sr=(1-h)*SSBo/(5*h-1);// Reclutamiento


// -----------------primer a�o
// genero una estructura inicial en torno a Z del primer a�o;
  Rnpre = mfexp(log_Rmed_pre);

  Neq=pre*Rnpre;
  for (j=1;j<=nper;j++)
   {
   Neq = (Neq*exp(-1.*M))*T + pre*exp(log_Rmed_pre+log_desv_No(j));
   }
  N(1)=Neq;
  NMD(1)=elem_prod(elem_prod(N(1),mfexp(-dt(1)*Z(1))),msex);
  BD(1)=sum(elem_prod(Wmed,NMD(1)));

// --------------------dinamica anual
  Rpop = mfexp(log_Rmed);
  Rpred(1) = Reclutas(1);

  for (i=2;i<=ntime;i++){

  Rpred(i)=Reclutas(i);

  if(i>minedad){

  Rpred(i)=(alfa_sr*BD(i-minedad)/(beta_sr+BD(i-minedad)));
  Reclutas(i)=Rpred(i)*mfexp(log_desv_Rt(i)); }

  N(i)=(elem_prod(N(i-1),S(i-1)))*T+pre*Reclutas(i);
  NMD(i)=elem_prod(elem_prod(N(i),mfexp(-dt(1)*Z(i))),msex);
  BD(i)=sum(elem_prod(Wmed,NMD(i)));

  } //


FUNCTION Eval_biomasas

 NMD=elem_prod(N,mfexp(-dt(1)*Z));
 NMD=elem_prod(NMD,outer_prod(Unos_year,msex));
 NVflo=elem_prod(elem_prod(N,mfexp(-dt(2)*(Z))),Sel);
 NVcru=elem_prod(elem_prod(N,mfexp(-dt(3)*(Z))),Selc);

// vectores de biomasas derivadas
// BD=Wmed*NMD;
 BMflo=Wmed*trans(NVflo);
 BMcru=Wmed*trans(NVcru);
 BMmph=Wmed*trans(NVmph);
 BT=Wmed*trans(N);

 BDo=sum(elem_prod(No,Wmed));
 RPRlp=BD/SSBo;


FUNCTION Eval_capturas_predichas

// matrices de capturas predichas por edad y año
 pred_Ctot=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));

// vectores de desembarques predichos por a�o
 pred_Desemb=Wmed*trans(pred_Ctot);

// matrices de proporcion de capturas por talla y a�o
 pobs=elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_tallas));
 ppred=elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot+1e-10),Unos_tallas));

 pobsc=elem_div(Ncru,outer_prod(rowsum(Ncru+1e-10),Unos_tallas));
 ppredc=elem_div(NVcru,outer_prod(rowsum(NVcru+1e-10),Unos_tallas));

 Lmed_pred=Tallas*trans(ppred);
 Lmed_obs=Tallas*trans(pobs);

 Lmed_predc=Tallas*trans(ppredc);
 Lmed_obsc=Tallas*trans(pobsc);


FUNCTION Eval_indices

   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques;j++){
              if (yrs(i)>=yqbloques(j)){
                 pred_CPUEI(i)=exp(log_qfloI(j))*pow(BMflo(i),exp(log_b));}
       }
   }

  for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques;j++){
              if (yrs(i)>=yqbloques(j)){
                 pred_CPUEA(i)=exp(log_qfloA(j))*pow(BMflo(i),exp(log_b));}
       }
   }

   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloquesc;j++){
              if (yrs(i)>=yqbloquesc(j)){
                 pred_Bcru(i)=exp(log_qcru(j))*BMcru(i);}
       }
   }


  for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloquesmph;j++){
              if (yrs(i)>=yqbloquesmph(j)){
                 pred_mph(i)=exp(log_qmph(j))*BD(i);}
       }
   }

  // INDICADORES DE REDUCCI�N DEL STOCK
  //RPRdin  = elem_div(BD,BDo);                       // RPR BDspr_t, din�mico
  //RPRequ  = BD/Bspro;                               // RPR con BDspro
  // RPRequ2 = BD/Bo;                                 // RPR con Bo proxy
  RPRequ3 = BD/(SSBo*0.55);                            // Raz�n para diagrama de fase
  Frpr    = exp(log_F)/0.47;

FUNCTION Eval_logverosim
// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
 int i;

 suma1=0; suma2=0; suma3=0; suma4=0; penalty=0;

 for (i=1;i<=ntime;i++)
 {
  if (CPUEI(i)>0){
    suma1+=square(log(CPUEI(i)/pred_CPUEI(i))*1/cv_index(2,i));}

  if (CPUEA(i)>0){
    suma2+=square(log(CPUEA(i)/pred_CPUEA(i))*1/cv_index(3,i));}

  if (Bcru(i)>0){
    suma3+=square(log(Bcru(i)/pred_Bcru(i))*1/cv_index(4,i));}

   if (mph(i)>0){
    suma4+=square(log(mph(i)/pred_mph(i))*1/cv_index(5,i));}
 }

FUNCTION Eval_funcion_objetivo

 suma5=0; suma6=0; penalty=0;

 likeval(1)=0.5*suma1;//CPUEI
 likeval(2)=0.5*suma2;//CPUEA
 likeval(3)=0.5*suma3;//Bcru
 likeval(4)=0.5*suma4;//mph

 likeval(5)=0.5*norm2(elem_div(log(elem_div(Desemb,pred_Desemb)),cv_index(1)));// desemb

 for (int i=1;i<=ntime;i++){
 suma5+=-nm(i)*sum(elem_prod(pobs(i),log(ppred(i)+1e-10)));
 suma6+=-nmc(i)*sum(elem_prod(pobsc(i),log(ppredc(i)+1e-10)));
 }

 likeval(6)=suma5;//
 likeval(7)=suma6;//

// lognormal Ninicial y Reclutas
 if(active(log_desv_Rt)){
 likeval(8)=1./(2*square(sigmaR))*norm2(log_desv_Rt);}

 if(active(log_desv_No)){
 likeval(9)=1./(2*square(sigmaR))*norm2(log_desv_No);}

 if(active(log_Lr)){
 likeval(10)=1./(2*square(cv_Par_bio(3)))*square(log_Lr-log_Lr_prior);}

  if (active(log_F)){
  pF=1000*norm2(log_F-mean(log_F));}

  penalty(1)=0.5/square(cv_q(1))*norm2(log_qfloI-log_prior_qfI);
  penalty(2)=0.5/square(cv_q(2))*norm2(log_qfloA-log_prior_qfA);
  penalty(3)=0.5/square(cv_q(3))*norm2(log_qcru-log_prior_qc);
  penalty(4)=0.5/square(cv_q(4))*norm2(log_qmph-log_prior_qmph);

  f=festim*(sum(likeval)+sum(penalty)+pF);

  if(last_phase){
  f=festim*(sum(likeval)+sum(penalty));}



FUNCTION  Eval_CTP

//-----------------------------------------------------------------


  for (int i=1;i<=npbr;i++){ // ciclo de PBR

  Np=N(ntime);
  Sp=S(ntime);


  //Fpbr=F(ntime)*pbr(i);// usa multiplicador pbr
  Fpbr=Sel(ntime)*pbr(i);// este usa la selectividad x Fpbr

  Zpbr=Fpbr+M;

  CTPp = elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Np));
  Ypact(i) = sum(elem_prod(CTPp,Wmed));
  for (int j=1;j<=ntime_sim;j++){ // ciclo de a�os


  if(j<=minedad){
   // Np=(elem_prod(Np,Sp))*T+pre*(alfa_sr*BD(ntime-minedad+1)/(beta_sr+BD(ntime-minedad+1)));} //Estima CTP con R_last
   //Np=(elem_prod(Np,Sp))*T+pre*(mfexp(log_Rmed))*2;} // Estima CTP con R_med*2 (R alto)
   Np=(elem_prod(Np,Sp))*T+pre*(mfexp(log_Rmed));} // Estima CTP con R_med

  if(j>minedad){
  Np=(elem_prod(Np,Sp))*T+pre*(alfa_sr*Bp(i,j-minedad)/(beta_sr+Bp(i,j-minedad)));} //

  //Rpp(i,j)=(alfa_sr*BD(ntime-minedad+1)/(beta_sr+BD(ntime-minedad+1)));

  Bp(i,j)=sum(elem_prod(elem_prod(Np,exp(-dt(1)*Zpbr)),elem_prod(msex,Wmed)));
  CTPp=elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Np));
  Yp(i,j)=sum(elem_prod(CTPp,Wmed));
  Sp=exp(-1.*Zpbr);
  }}


REPORT_SECTION

 report << "Years" << endl;
 report << yrs << endl;
 report << "Bcru_obs" << endl;
 report << Bcru << endl;
 report << "Bcru_pred" << endl;
 report << pred_Bcru << endl;
 report << "cpue1_obs" << endl;
 report << CPUEI << endl;
 report << "pred_cpue1" << endl;
 report << pred_CPUEI << endl;
 report << "cpue2_obs" << endl;
 report << CPUEA << endl;
 report << "pred_cpue2" << endl;
 report << pred_CPUEA << endl;
 report << "MPH_obs" << endl;
 report << mph << endl;
 report << "MPH_pred" << endl;
 report << pred_mph << endl;
 report << "Desemb_obs" << endl;
 report << Desemb << endl;
 report << "Desemb_pred" << endl;
 report << pred_Desemb << endl;
 report << "Lmf_obs" << endl;
 report << Lmed_obs << endl;
 report << "Lmf_pred" << endl;
 report << Lmed_pred << endl;
 report << "Lmc_obs" << endl;
 report << Lmed_obsc << endl;
 report << "Lmc_pred" << endl;
 report << Lmed_predc << endl;
 report << "Biomasa_desovante" << endl;
 report << BD << endl;
 report << "Biomasa_total" << endl;
 report << BT << endl;
 report << "Biomasa_explotable" << endl;
 report << BMflo << endl;
 report << "N_explotable" << endl;
 report << pred_Ctot << endl;
 report << "Reclutamiento" << endl;
 report << Reclutas<< endl;
 report << "Rpred" << endl;
 report << Rpred<< endl;
 report << "F" << endl;
 report << exp(log_F) << endl;
 report<<"Tallas"<<endl;
 report<<Tallas<<endl;
 report<<"Wmed"<<endl;
 report<<Wmed<<endl;
 report<<"Msex"<<endl;
 report<<msex<<endl;
 report<<"Abundancia_talla"<<endl;
 report<<N<<endl;
 report<<"Selflo_talla"<<endl;
 report<<Sel<<endl;
 report<<"Selcru_talla"<<endl;
 report<<Selc<<endl;
 report << "Propfl_obs" << endl;
 report << pobs<< endl;
 report << "Propfl_pred" << endl;
 report << ppred<< endl;
 report << "Propcru_obs" << endl;
 report << pobsc<< endl;
 report << "Propcru_pred" << endl;
 report << ppredc<< endl;
 report << "BD_virgen_anual" << endl;
 report << BDo << endl;
 report << "BD_virgen_LP" << endl;
 report << SSBo << endl;
 report << "Reduccion_LP " << endl;
 report << RPRlp << endl;
 report << "Talla_media_por_grupo" << endl;
 report << Lesp << endl;
 report <<  "desvest_por_grupo" << endl;
 report << sigmaL << endl;
 report << "Fun_rec_talla" << endl;
 report << pre<< endl;
 report << "MatrizTrans" << endl;
 report << T << endl;
 report << "bCPUE  Lr  Sr  beta  h " << endl;
 report << exp(log_b)<<" "<<exp(log_Lr)<<" "<<exp(log_sr)<<" "<<exp(log_beta)<<" "<<h<< endl;


//-------------------------------------------------------------------
// ESTIMA nm y CV

  suma1=0; suma2=0; nm1=1; cuenta1=0;

  for (int i=1;i<=ntime;i++){ //

   if (sum(pobs(i))>0){
      suma1=sum(elem_prod(ppred(i),1-ppred(i)));
      suma2=norm2(pobs(i)-ppred(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

 report << "nm_flota_cru" <<endl;
 report <<pow(nm1,1/cuenta1)<< endl;


 suma1=0; suma2=0;nm1=1;cuenta1=0;

  for (int i=1;i<=ntime;i++){ //

   if (sum(pobs(i))>0){
      suma1=sum(elem_prod(ppredc(i),1-ppredc(i)));
      suma2=norm2(pobsc(i)-ppredc(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

 report <<pow(nm1,1/cuenta1)<< endl;

 report << "BD_proy" << endl; //biomasa desovante proyectada Flast"
 report << Bp << endl;
 report << "Capt_proy" << endl; // Capturas proyectadas para cada Fpbr
 report << Yp << endl;
 report << "Captura_act" << endl;  //Captura proyectadas a�o en curso
 report << Ypact << endl;
 report << "LIKE" << endl;
 report << likeval << endl;
 report << "Prioris" << endl;
 report << penalty << endl;
 report << "log_Rmed" << endl;
 report << priorlog_Rmed << endl;
 report << "Np" << endl;
 report << Np << endl;
 report << "PPp" << endl;
 report << pre*(mfexp(log_Rmed)) << endl;
 report << "Nvp" << endl;
 report << CTPp << endl;


FINAL_SECTION

 time(&finish);
 elapsed_time=difftime(finish,start);
 hour=long(elapsed_time)/3600;
 minute=long(elapsed_time)%3600/60;
 second=(long(elapsed_time)%3600)%60;
 cout<<endl<<endl<<"*********************************************"<<endl;
 cout<<"--Start time:  "<<ctime(&start)<<endl;
 cout<<"--Finish time: "<<ctime(&finish)<<endl;
 cout<<"--Runtime: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
 cout<<"*********************************************"<<endl;
