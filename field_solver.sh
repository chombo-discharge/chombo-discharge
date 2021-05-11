# Find and replace everything field_solver_multigrid.bc_x_lo by field_solver_multigrid.bc.x.lo
for i in `find . -name "*.options" -type f`; do
    sed -i 's/field_solver_multigrid.bc_x_low/field_solver_multigrid.bc.x.lo /g' $i
    sed -i 's/field_solver_multigrid.bc_x_high/field_solver_multigrid.bc.x.hi  /g' $i

    sed -i 's/field_solver_multigrid.bc_y_low/field_solver_multigrid.bc.y.lo /g' $i
    sed -i 's/field_solver_multigrid.bc_y_high/field_solver_multigrid.bc.y.hi  /g' $i

    sed -i 's/field_solver_multigrid.bc_z_low/field_solver_multigrid.bc.z.lo /g' $i
    sed -i 's/field_solver_multigrid.bc_z_high/field_solver_multigrid.bc.z.hi  /g' $i


    sed -i 's/= neumann    /= neumann 0.0/g' $i
    sed -i 's/= dirichlet_ground/= dirichlet 0.0   /g' $i
    sed -i 's/= dirichlet_live/= dirichlet 1.0 /g' $i

    # Delete old doc description. 
    sed -i 's/BC type. \"neumann\", \"dirichlet_ground\", \"dirichlet_live\", or \"robin\"/Bc type./g' $i
done

for i in `find . -name "*.inputs" -type f`; do
    sed -i 's/field_solver_multigrid.bc_x_low/field_solver_multigrid.bc.x.lo /g' $i
    sed -i 's/field_solver_multigrid.bc_x_high/field_solver_multigrid.bc.x.hi  /g' $i

    sed -i 's/field_solver_multigrid.bc_y_low/field_solver_multigrid.bc.y.lo /g' $i
    sed -i 's/field_solver_multigrid.bc_y_high/field_solver_multigrid.bc.y.hi  /g' $i

    sed -i 's/field_solver_multigrid.bc_z_low/field_solver_multigrid.bc.z.lo /g' $i
    sed -i 's/field_solver_multigrid.bc_z_high/field_solver_multigrid.bc.z.hi  /g' $i


    sed -i 's/= neumann /= neumann 0.0/g' $i
    sed -i 's/= dirichlet_ground/= dirichlet 0.0   /g' $i
    sed -i 's/= dirichlet_live/= dirichlet 1.0 /g' $i

    # Delete old doc description. 
    sed -i 's/BC type. \"neumann\", \"dirichlet_ground\", \"dirichlet_live\", or \"robin\"/Bc type./g' $i
done
