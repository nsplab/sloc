#!/bin/bash

data=${HOME}/meshes

#air_scalp_layer=1
#scalp_skull_layer=2
#skull_csf_layer=3
#csf_brain_layer=4
#brain_csf_layer=5
#brain_plasma_layer=6

air_scalp_layer_id=1
scalp_skull_layer_id=2
skull_brain_layer_id=3
brain_csf_layer_id=4
brain_plasma_layer_id=5

air_scalp_layer=air/scalp
scalp_skull_layer=scalp/skull
skull_brain_layer=skull/brain
brain_csf_layer=brain/csf
brain_plasma_layer=brain/plasma

scalp_mat=air_scalp_layer_id
skull_mat=scalp_skull_layer_id
brain_mat=skull_brain_layer_id
ventricles_mat=brain_csf_layer_id
artery_mat=brain_plasma_layer_id

scalp_stl=${data}/scalp11.stl
skull_stl=${data}/skull11.stl
brain_stl=${data}/brain11.stl
ventricles_stl=${data}/ventricles11.stl
artery_stl=${data}/artery11.stl

scalp_ucd=./scalp.ucd
skull_ucd=./skull.ucd
brain_ucd=./brain.ucd
ventricles_ucd=./ventricles.ucd
artery_ucd=./artery.ucd

if [ ! -f ./head.surf.ucd ]; then

    for var in scalp skull brain ventricles artery ; do
        var_stl=${var}_stl
        var_ucd=${var}_ucd
        var_mat=${var}_mat
        mat_id=${!var_mat}
        #stl_to_ucd -i ${!var_stl} -o ${var}.tri.ucd
        #tri_to_quad -i $f.tri.ucd -o ${var}.quad.ucd
        #set_material_id -m ${!mat_id} -i ${var}.quad.ucd -o ${!var_ucd}
        echo ----
        stl_to_ucd ${!var_stl} ${var}.tri.ucd
        echo ----
        tri_to_quad ${var}.tri.ucd ${var}.quad.ucd
        echo ----
        set_material_id ${!mat_id} ${var}.quad.ucd ${!var_ucd}
    done

    #merge_ucd_files ${scalp_ucd} ${skull_ucd} ${brain_ucd} ${ventricles_ucd} ${artery_ucd} -o head.surf.ucd 

    # do a pairwise merge for now
    echo ----
    merge_ucd_files ${scalp_ucd} ${skull_ucd} foo.ucd
    echo ----
    merge_ucd_files foo.ucd ${brain_ucd} bar.ucd
    echo ----
    merge_ucd_files bar.ucd ${ventricles_ucd} foo.ucd
    echo ----
    merge_ucd_files foo.ucd ${artery_ucd} bar.ucd
    echo ----
    rm -f foo.ucd
    mv bar.ucd head.surf.ucd

    # clean up temp files
    rm -f *.tri.ucd *.quad.ucd
fi

#create_material_layers.py \
#    material_air=3e-15 \
#    material_scalp=0.275 \
#    material_skull=0.0132 \
#    material_csf=1.79 \
#    material_brain=0.40 \
#    material_plasma=0.667 \
#    layer_1=csf/brain \
#    layer_2=skull/brain \
#    layer_3=scalp/skull \
#    layer_4=air/scalp \
#    output=./head.sigma

create_material_layers.py \
    material_air=3e-15 \
    material_scalp=0.275 \
    material_skull=0.0132 \
    material_csf=1.79 \
    material_brain=0.40 \
    material_plasma=0.667 \
    layer_1=${air_scalp_layer} \
    layer_2=${scalp_skull_layer} \
    layer_3=${skull_brain_layer} \
    layer_4=${brain_csf_layer} \
    layer_5=${brain_plasma_layer} \
    output=./head.sigma

