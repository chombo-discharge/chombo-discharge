Using VisIt on clusters {#visualization-cluster}
=======================

To use VisIt's parallel computing capabilities to its fullest, its best to run it in server mode. To do this, _you should have matching versions of VisIt both locally and remotely_.

To set up remote visualization, you will need

1. The host address.
2. The path to VisIt on the server
3. Optionally have to modify your .bashrc (or equivalent) so that VisIt is loaded on ssh login.

Below, we show how to set up a host profile for fram.sigma2.no

Remote visualization on fram.sigma2.no
--------------------------------------

We will use the remote 2.13.0 remote installation of VisIt on fram. Firstly, modify your .bashrc file so that your login automatically loads the necessary modules. On fram, this is done by appending the following line to your .bashrc file

      module load VisIt/2.13.0-intel-2017a

Next, we will set up the host profile on your local installation. Run VisIt2.13.0 locally

1. Go to Options->Host profiles
2. Create a new host profile
3. In the field 'Host name', use any host name. 
4. In the field 'Remote host name', type 'fram.sigma2.no'
5. In the field 'Path to VisIt installation', type '/cluster/software/VisIt/2.13.0-intel-2017a'
6. Check the 'Tunnel data connections through SSH' flag.
7. Navigate to the 'Launch profiles' tab.
8. Create a new profile and navigate to the 'Parallel' tab.
9. In the 'Parallel' tab, check the 'Launch parallel engine' box.

To open a file remotely. Select 'Open' and choose your new host under the 'Host' dropdown menu. You will be prompted for a password, after which you can navigate to your remote file. 