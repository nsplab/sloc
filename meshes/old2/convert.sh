#!/bin/bash -x

# -----------------------------------------------------------------------------

[ ! -f Bone3.ucd ] \
    && stl_to_ucd Bone3.stl Bone3-tri.ucd \
    && tri_to_quad Bone3-tri.ucd Bone3.ucd

[ ! -f Brain3.ucd ] \
    && stl_to_ucd Brain3.stl Brain3-tri.ucd \
    && tri_to_quad Brain3-tri.ucd Brain3.ucd

[ ! -f Artery3.ucd ] \
    && stl_to_ucd Artery3.stl Artery3-tri.ucd \
    && tri_to_quad Artery3-tri.ucd Artery3.ucd

[ ! -f Vessels3.ucd ] \
    && stl_to_ucd Vessels3.stl Vessels3-tri.ucd \
    && tri_to_quad Vessels3-tri.ucd Vessels3.ucd

# -----------------------------------------------------------------------------

[ ! -f Bone2.ucd ] \
    && stl_to_ucd Bone2.stl Bone2-tri.ucd \
    && tri_to_quad Bone2-tri.ucd Bone2.ucd

[ ! -f Brain2.ucd ] \
    && stl_to_ucd Brain2.stl Brain2-tri.ucd \
    && tri_to_quad Brain2-tri.ucd Brain2.ucd

[ ! -f Artery2.ucd ] \
    && stl_to_ucd Artery2.stl Artery2-tri.ucd \
    && tri_to_quad Artery2-tri.ucd Artery2.ucd

[ ! -f Vessels2.ucd ] \
    && stl_to_ucd Vessels2.stl Vessels2-tri.ucd \
    && tri_to_quad Vessels2-tri.ucd Vessels2.ucd

# -----------------------------------------------------------------------------

[ ! -f Bone.ucd ] \
    && stl_to_ucd Bone.stl Bone-tri.ucd \
    && tri_to_quad Bone-tri.ucd Bone.ucd

[ ! -f Brain.ucd ] \
    && stl_to_ucd Brain.stl Brain-tri.ucd \
    && tri_to_quad Brain-tri.ucd Brain.ucd

[ ! -f Brain1a.ucd ] \
    && stl_to_ucd Brain1a.stl Brain1a-tri.ucd \
    && tri_to_quad Brain1a-tri.ucd Brain1a.ucd

[ ! -f Artery.ucd ] \
    && stl_to_ucd Artery.stl Artery-tri.ucd \
    && tri_to_quad Artery-tri.ucd Artery.ucd

[ ! -f Vessels.ucd ] \
    && stl_to_ucd Vessels.stl Vessels-tri.ucd \
    && tri_to_quad Vessels-tri.ucd Vessels.ucd

# -----------------------------------------------------------------------------
