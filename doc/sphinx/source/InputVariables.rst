.. _Chap:InputVariables:

Input variables
---------------

Generally, the coding style for input variables is to use the class name as a prefix (where :ref:`Chap:amr_mesh` is an exception) and the variable as a suffix. All letters are lower-case. For example::

   plasma_engine.max_steps = 10

To pass input variables into PlasmaC, we generally refrain from hard-coding. Instead, we use Chombo's ParmParse class, which is used in the following way:

.. code-block:: c++

   Real my_variable;
   ParmParse pp("prefix");
   pp.get("suffix", my_variable);

The above code segment will attept to fetch an input line **prefix.suffix** and place it in *my_variable*. Note that the specification of **prefix.suffix** should be of the same type as *my_variable* (float in this case). For this example, passing

.. code-block:: bash

		mpirun -np 32 <my_application> <my_input_file> prefix.suffix = foo

will throw an error. Please refer to :ref:`Chap:NewSimulations` for more details. 
