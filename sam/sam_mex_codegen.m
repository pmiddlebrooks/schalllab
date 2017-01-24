function sam_mex_codegen(task)
% Converts some SAM functions to MEX files
%  
% SYNTAX 
% SAM_MEX_CODEGEN('sam_sim_trial'); 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Thu 25 Jul 2013 15:29:37 CDT by bram 
% $Modified: Tue 24 Sep 2013 11:32:15 CDT by bram

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. SPECIFY GENERAL SETTINGS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.MATLABSourceComments = true;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. SPECIFY WHAT FUNCTION TO DO
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

switch lower(task)
  case 'sam_sim_trial'
    codegen_sam_sim_trial;
  case 'sam_spec_timing_diagram'
    codegen_sam_spec_timing_diagram
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3. NESTED FUNCTIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function codegen_sam_spec_timing_diagram

  % This function compiles the sam_spec_timing_diagram.m file to a MEX file
  
  ons           = coder.typeof(0,[5e1,5e1],[1,1]);
  dur           = coder.typeof(0,[5e1,5e1],[1,1]);
  v             = coder.typeof(0,[5e1,1],[1,0]);
  eta           = coder.typeof(0,[5e1,1],[1,0]);
  se            = coder.typeof(0,[5e1,5e1],[1,1]);
  dt            = coder.typeof(0,[1,1],[0,0]);
  tWindow       = coder.typeof(0,[1,2],[0,0]);
  
  codegen sam_spec_timing_diagram -config cfg -args {ons,dur,v,eta,se,dt,tWindow}
  
end

function codegen_sam_sim_trial

  % This function compiles the sam_sim_trial.m file to a MEX file
  u           = coder.typeof(0,[5e1,1e6],[1,1]);
  A           = coder.typeof(0,[5e1,5e1],[1,1]);
  B           = coder.typeof(0,[5e1,5e1,1e3],[1,1,1]);
  C           = coder.typeof(0,[5e1,1e3],[1,1]);
  D           = coder.typeof(0,[5e1,5e1,5e1],[1,1,1]);
  Sin         = coder.typeof(0,[5e1,5e1],[1,1]);
  Z0          = coder.typeof(0,[5e1,1],[1,0]);
  ZC          = coder.typeof(0,[5e1,1],[1,0]);
  ZLB         = coder.typeof(0,[5e1,1],[1,0]);
  dt          = coder.typeof(0,[1,1],[0,0]);
  tau         = coder.typeof(0,[1,1],[0,0]);
  T           = coder.typeof(0,[1,1e6],[0,1]);
  terminate   = coder.typeof(true,[5e1,1],[1,0]);
  blockInput  = coder.typeof(true,[5e1,5e1],[1,1]);
  latInhib    = coder.typeof(true,[5e1,5e1],[1,1]);
  n           = coder.typeof(0,[1,1],[0,0]);
  m           = coder.typeof(0,[1,1],[0,0]);
  p           = coder.typeof(0,[1,1],[0,0]);
  t           = coder.typeof(0,[1,1],[0,0]);
  rt          = coder.typeof(0,[5e1,1],[1,0]);
  resp        = coder.typeof(true,[5e1,1],[1,0]);
  z           = coder.typeof(0,[5e1,1e6],[1,1]);

  codegen sam_sim_trial -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_crace_irace_nomodbd -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_crace_ibi_nomodbd -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_crace_ili_nomodbd -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cffi_irace_nomodbd -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cffi_ibi_nomodbd -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cffi_ili_nomodbd -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cli_irace_nomodbd -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cli_ibi_nomodbd -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cli_ili_nomodbd -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}

  codegen sam_sim_trial_inpdepnoise -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_crace_irace_nomodbd_inpdepnoise -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_crace_ibi_nomodbd_inpdepnoise -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_crace_ili_nomodbd_inpdepnoise -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cffi_irace_nomodbd_inpdepnoise -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cffi_ibi_nomodbd_inpdepnoise -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cffi_ili_nomodbd_inpdepnoise -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cli_irace_nomodbd_inpdepnoise -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cli_ibi_nomodbd_inpdepnoise -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}
  codegen sam_sim_trial_cli_ili_nomodbd_inpdepnoise -config cfg -args {u,A,B,C,D,Sin,Z0,ZC,ZLB,dt,tau,T,terminate,blockInput,latInhib,n,m,p,t,rt,resp,z}

end

end