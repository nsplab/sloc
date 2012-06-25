#!/bin/bash


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

if [ ! -f ./head.surf.ucd ]; then

    #for f in scalp skull brain artery vacuoles; do
    for f in scalp skull brain; do
        file=${f}_file
        layer=${!file}
        mat_id=${!layer}
        stl_to_ucd -i $f.stl -o $f.tri.ucd
        tri_to_quad -i $f.tri.ucd -o $f.quad.ucd
        set_material_id -m ${mat_id} -i $f.quad.ucd -o $f.ucd
    done

    #merge_ucd_files head.surf.ucd scalp.ucd skull.ucd brain.ucd artery.ucd vacuoles.ucd
    merge_ucd_files head.surf.ucd scalp.ucd skull.ucd brain.ucd

    #rm -f *.tri.ucd *.quad.ucd
fi

create_material_layers.py \
    material_air=3e-15 \
    material_scalp=0.275 \
    material_skull=0.0132 \
    material_csf=1.79 \
    material_brain=0.40 \
    material_plasma=0.667 \
    layer_1=csf/brain \
    layer_2=skull/brain \
    layer_3=scalp/skull \
    layer_4=air/scalp \
    output=head.sigma

