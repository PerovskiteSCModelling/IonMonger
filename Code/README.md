The scripts are divided into three folders:

- Common contains the scripts to create the necessary inputs, run the solver and process the raw solution
- Plotting contains scripts to plot a saved solution
- Solver contains the functions called by FE_solve.m to obtain a solution using a method based on finite elements

The Tests/ folder is for developers and contains a set of tests to check the consistency of future updates. The tests are written in the style of MATLAB's built-in testing framework. To run the tests use the following commands:
```
reset_path;
runtests('Code\Tests\test_IonMonger.m');
```
