% This script provides an examples of how to use the SAFE Toolbox in
% combination with a model that does not run under matlab.
%
% STEPS
%
% 1) In Matlab: use SAFE sampling functions to sample the input space. 
% Then save the input samples in a text file that will be passed on to the
% external model.
%
% 2) Outside Matlab [not shown in this workflow]: run the model against
% each input sample and save the corresponding output into a text file.
%
% 3) In Matlab: load the text file and use SAFE post-processing functions
% to compute sensitivity indices.

% This script prepared by Francesca Pianosi and Fanny Sarrazin
% University of Bristol, 2014
% mail to: francesca.pianosi@bristol.ac.uk

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%% Step 1

% a) Define the input feasible space

M = 5 ; % Number of inputs
DistrFun  = 'unif'  ; % Distribution (same for all inputs)
DistrPar  = { [ 0 400 ]; ...
              [ 0   2 ]; ...
              [ 0   1 ]; ...
              [ 0 0.1 ]; ...
              [ 0.1 1 ] } ; % Ranges

% b) Perform sampling

% For instance, One-At-the-Time sampling (see workflow about EET to learn
% more about available options for OAT):
r = 30 ; % number of sampling points
SampStrategy = 'lhs' ; % Latin Hypercube
design_type = 'radial'; % Type of design
X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);

% Or, All-At-the-Time sampling (again, see workflow about RSA for more
% options):
SampStrategy = 'lhs' ;
N = 3000 ; % Number of samples
X = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,N);

% c) Save to file:

% Choose name for the input file:
file_name = 'Input_samples.txt' ;

% Set current directory to the directory where the input file will be saved:
cd([ my_dir '/example/external'])

% Write to file:
save(file_name,'X','-ascii') ; % With this command, data will be saved in 
% exponential notation with 8 digits.

% If you don't like that format, you can use the fprintf command, e.g.:
formattype = '%3.0f \t %1.2f \t %1.2f \t %1.3f \t %1.3f \n' ;
fid = fopen(file_name,'w')  ;
fprintf(fid,formattype,X');
fclose(fid);
% Or, let the function choose the 'more compact' format (but in this case
% double-check that numbers in the file have the required precision):
formattype = '' ; for i=1:M-1; formattype = [ formattype '%g\t' ]; end; formattype = [ formattype '%g\n' ];
fid = fopen(file_name,'w')  ;
fprintf(fid,formattype,X');
fclose(fid);

% 'fprintf' is also useful to add an header to the file:
header = 'This file created by XXX' ;
formattype = '%3.0f \t %1.2f \t %1.2f \t %1.3f \t %1.3f \n' ;
fid = fopen(file_name,'w')  ;
fprintf(fid,'%s\n',header);
fprintf(fid,formattype,X');
fclose(fid);

%% Step 2

% Run the model outside matlab and generate the output file


%% Step 3

% Set current directory to the directory where the output file was saved:
cd([ my_dir '/example/external'])

% Output file name:
output_file = 'Output_samples.txt' ;

% Load data from file:
Y = load(output_file,'-ascii') ;

% From now on, just use all the functions in RSA, EET, VBSA, etc.



