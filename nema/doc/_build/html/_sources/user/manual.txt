Launching Nema
################

How to run Nema from the test directory
----------------------------------------

Nema needs several input arguments:
    * *inputdir*: input directory 
    * *outputdir*: output directory where nema write the results
    * *outf*: not used variable


When you run **nema**, be sure that you have created the following directories in the input directory:
    * parameters
    * DrivingVariables
    * State0dd

Please, create also an outputs directory in the outputdir.

.. code-block:: bash

    my@alinea:~nema/test$ ls
    DrivingVariables  outputs  parameters  State0dd

    # run nema program in test directory
    my@alinea:~nema/test$ nema

