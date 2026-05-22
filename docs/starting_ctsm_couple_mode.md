## To get started for ctsm coupling mode

User interested in using mizuRoute with CTSM is referred to CESM's user guide. Here, quick guide is provided.  

1. Obtain CTSM code from [github](https://github.com/ESCOMP/CTSM/tree/master)

2. Create the case

   ```bash
   cd cime/scripts
   ```

   ```bash
   ./create_newcase --case <testcase> --mach derecho --res f09_f09_rHDMAlk_mg17 -compset I2000Clm60SpMizGs 
   ```
   (`./create_newcase -help --` to get help on the script)

3. Setup the case

   ```bash
   cd <testcase>
   ```

   ```bash
   ./xmlchange id1=val1,id2=val2  # to make changes to any settings in the env_*.xml files
   ./case.setup
   ```
   (./case.setup -help -- to get help on the script)

3. Add any namelist changes to the `user_nl_*` files

   ```bash
   $EDITOR user_nl_*
   ```

4. Compile the code

   ```bash
   ./case.build
   ```

5. Submit the run

   ```bash
   ./case.submit
   ```


