# GOMC - GPU Optimized Monte Carlo

## Building GOMC on GNU/Linux, macOS, or Cygwin:

  1. Clone or download our code from GitHub:
      ```bash
      git clone https://github.com/GOMC-WSU/GOMC.git
      ```
  2. Go into the GOMC directory: 
      ```bash
      cd GOMC
      ```
  3. Checkout to `SDSC` branch.
      ```bash
      git checkout SDSC
      ```
  5. Give execution permission: 
      ```bash
      chmod u+x metamake.sh
      ```
  6. Run metamake file:
      ```bash
      ./metamake.sh
      ```
  7. Last step should generate all the executables in `bin` directory

  You can either copy the executable into your example directory or add `bin` directory to your path. We recommend the latter.

## Testing
  An example has been included in the `example` directory. You can run using the following command:
  ```bash
  GOMC_GPU_NPT in.conf | tee gpu.log
  ```

  Using `tee` you should be able to see your output and redirect it to the file `gpu.log` as well. To confirm correctness you can `diff` with our reference output.
  ```bash
  diff gpu.log ref.log
  ```

  If you do not see a bunch of scientific numbers you are good to go!
