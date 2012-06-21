#!/bin/bash

# -----------------------------------------------------------------------------

if [ ! -f ./head.surf.ucd ]; then

    #air_scalp_layer=1
    #scalp_skull_layer=2
    #skull_csf_layer=3
    #csf_brain_layer=4
    #brain_csf_layer=5
    #brain_plasma_layer=6

    air_scalp_layer=1
    scalp_skull_layer=2
    skull_brain_layer=3
    brain_plasma_layer=4
    brain_csf_layer=5

    scalp_file=air_scalp_layer
    skull_file=scalp_skull_layer
    brain_file=skull_brain_layer
    artery_file=brain_plasma_layer
    vacuoles_file=brain_csf_layer

    for f in scalp skull brain artery vacuoles; do
        file=${f}_file
        layer=${!file}
        mat_id=${!layer}
        stl_to_ucd -i $f.stl -o $f.tri.ucd
        tri_to_quad -i $f.tri.ucd -o $f.quad.ucd
        set_material_id -m ${mat_id} -i $f.quad.ucd -o $f.ucd
    done

    merge_ucd_files \
        -o head.surf.ucd \
        scalp.ucd skull.ucd brain.ucd artery.ucd vacuoles.ucd
fi

create_material_layers
    --material-air=3e-15 \
    --material-scalp=0.275 \
    --material-skull=0.0132 \
    --material-csf=1.79 \
    --material-brain=0.40 \
    --material-plasma=0.667 \
    --layer-1=csf/brain \
    --layer-2=skull/brain \
    --layer-3=scalp/skull \
    --layer-4=air/scalp \
    --output=head.sigma

# -----------------------------------------------------------------------------

create_dipoles \
    --output=head.dipoles \
    0,0,0 \
    1,1,1 \
    ...

select_electrodes \
    --input-mesh=head.surf.ucd \
    --electrodes-vtk=head.ec1.vtk \
    --electrodes=head.ec1 \
    1/30 4/35 3/30

create_forward_params \
    --debug=true \
    --verbose=true \
    --surface-mesh=head.surf.ucd \
    --material-data=head.sigma \
    --dipoles=head.dipoles \
    --output-vtk=head.vtk \
    --output-phi=head_phi.dat \
    head.fwd.prm

bem_forward_solver ... head.fwd.prm

measure_electrodes \
    --phi=./head_phi.dat \
    --electrodes=./head.ec1 \
    --snr=3.2 \
    --output-measurements=./head.em1

create_inverse_params \
    --debug=true \
    --verbose=false \
    --input-measurements=head.em1 \
    --output-sources=head.source.dipoles \
    --fwd-surface-mesh=head.surf.ucd \
    --fwd-material-data=head.sigma \
    --fwd-verbose=true \
    --search-verbose=true \
    --search-debug=true \
    --search-initial-point=0,0,0 \
    --search-initial-radius=10e-3 \
    --search-tolerance=1e-8 \
    --search-max-iterations=1000 \
    --search-reflection-coefficient=1.0 \
    --search-contraction-coefficient=0.5 \
    --search-expansion-coefficient=2.0 \
    --search-reduction-coefficient=0.5 \
    head.inv.prm

bem_inverse_solver ... head.inv.prm

