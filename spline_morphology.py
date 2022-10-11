
import bpy
import mathutils

#
# Object management utils.
#

def select_activate_only(objects=[]):
    for obj in bpy.data.objects:
        obj.select_set(False)
    bpy.context.view_layer.objects.active = None 
    for obj in objects:
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
    
def copy_object(base_object, collection_name=None):
    copy_object = base_object.copy()
    copy_object.data = base_object.data.copy()
    copy_object.animation_data_clear()
    if collection_name == None:
        bpy.context.collection.objects.link(copy_object)
    else:
        add_object_to_collection(copy_object, collection_name)
    return copy_object

# https://graphics.pixar.com/library/OrthonormalB/paper.pdf
# NOTE: n must be normalized!
def pixar_onb(_n):
    n = _n.normalized()
    t = mathutils.Vector((0,0,0))
    b = mathutils.Vector((0,0,0))
    if(n[2] < 0.0):
        a = 1.0 / (1.0 - n[2])
        b = n[0] * n[1] * a
        t = mathutils.Vector((1.0 - n[0] * n[0] * a, -b, n[0]))
        b = mathutils.Vector((b, n[1] * n[1] * a - 1.0, -n[1]))
    else:
        a = 1.0 / (1.0 + n[2])
        b = -n[0] * n[1] * a
        t = mathutils.Vector((1.0 - n[0] * n[0] * a, b, -n[0]))
        b = mathutils.Vector((b, 1 - n[1] * n[1] * a, -n[1]))
    return t, b

def transform_obj_rot_pos(base_obj, pos_vec=mathutils.Vector((0,0,0)), rot_vec=mathutils.Vector((0,0,0))):
    # Extract position and scale.
    #original_scale_vec = base_obj.matrix_world.to_scale()
    #original_mat_scale = mathutils.Matrix.Scale(0.5, 4, (0.0, 0.0, 1.0))
    #original_pos_vec = base_obj.matrix_world.to_translation()
    #original_mat_trans = mathutils.Matrix.Translation(original_pos_vec)
    # 
    # zero out curr rotation matrix first: https://blender.stackexchange.com/a/159992
    curr_rot_mat_inv = base_obj.matrix_basis.to_3x3().transposed().to_4x4()
    base_obj.matrix_basis = base_obj.matrix_basis @ curr_rot_mat_inv
    # Orient object using vector.
    new_rot_z = rot_vec.normalized()
    new_rot_x, new_rot_y = pixar_onb(new_rot_z)
    rot_basis = mathutils.Matrix((new_rot_x, new_rot_y, new_rot_z))
    rot_basis = rot_basis.transposed()
    rot_basis.resize_4x4()
    rot_mat = rot_basis.to_euler().to_matrix().to_4x4() # extract only rotation!
    base_obj.matrix_basis = base_obj.matrix_basis @ rot_mat # https://blender.stackexchange.com/questions/35125/what-is-matrix-basis
    # Transform object using vector.
    new_pos = pos_vec
    trans_mat = mathutils.Matrix.Translation(new_pos)
    base_obj.matrix_basis = trans_mat @ base_obj.matrix_basis

def create_instance(base_obj, pos_vec=mathutils.Vector((0,0,0)), rot_vec=mathutils.Vector((0,0,0)), collection_name=None):
    # Create instance.
    inst_obj = bpy.data.objects.new(base_obj.name+"_inst", base_obj.data)
    # Transform.
    transform_obj_rot_pos(inst_obj, pos_vec, rot_vec)
    # Store.
    if collection_name == None:
        bpy.context.collection.objects.link(inst_obj)
    else:
        create_collection_if_not_exists(collection_name)
        bpy.data.collections[collection_name].objects.link(inst_obj)
    return inst_obj

# https://blender.stackexchange.com/questions/220072/check-using-name-if-a-collection-exists-in-blend-is-linked-to-scene
def create_collection_if_not_exists(collection_name):
    if collection_name not in bpy.data.collections:
        new_collection = bpy.data.collections.new(collection_name)
        bpy.context.scene.collection.children.link(new_collection) #Creates a new collection

def add_object_to_collection(base_object, collection_name="collection"):
    create_collection_if_not_exists(collection_name)
    bpy.data.collections[collection_name].objects.link(base_object)

def move_object_to_collection(base_object, collection_name="collection"):
    create_collection_if_not_exists(collection_name)
    for curr_collection in base_object.users_collection:
        curr_collection.objects.unlink(base_object)
    bpy.data.collections[collection_name].objects.link(base_object)
    

#
# Math utils.
#

# Interpolate [a,b] using factor t.
def lerp(t, a, b):
    return (1.0 - t) * a + t * b

#
# Curve utils.
#

def remove_curve_points(curve_obj):
    # TODO
    pass

def subdivide_curve():
    # TODO
    pass

def create_bezier_curve_bops(name="bezier_curve", collection_name=None):
    select_activate_only([])
    bpy.ops.curve.primitive_bezier_curve_add(enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
    created_curve = bpy.context.selected_objects[0]
    created_curve.name = name
    if (collection_name != None):
        move_object_to_collection(created_curve, collection_name)
    return created_curve

def create_nurbs_path_curve_bops(name="bezier_curve", collection_name=None):
    select_activate_only([])
    bpy.ops.curve.primitive_nurbs_path_add(radius=1, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
    created_curve = bpy.context.selected_objects[0]
    created_curve.name = name
    if (collection_name != None):
        move_object_to_collection(created_curve, collection_name)
    return created_curve

# Perturb curve.
def perturb_curve(curve_obj, perturb_scale=1.0, perturb_strength=1.0, n_octaves=1, amplitude_scale=1.0, frequency_scale=1.0):
    curve_type = curve_obj.data.splines[0].type
    points = []
    if curve_type == "BEZIER":
        points = curve_obj.data.splines[0].bezier_points
    if curve_type == "POLY" or curve_type == "NURBS":
        points = curve_obj.data.splines[0].points
    for point in points:
        point_co = mathutils.Vector((point.co[0], point.co[1], point.co[2]))       
        trans_vec = mathutils.noise.turbulence_vector(
            point_co * perturb_scale, 
            n_octaves,
            False, #hard
            noise_basis='PERLIN_ORIGINAL',
            amplitude_scale=amplitude_scale,
            frequency_scale=frequency_scale) * perturb_strength
        new_point_co = point_co + trans_vec
        if curve_type == "BEZIER":
            point.co = (new_point_co[0], new_point_co[1], new_point_co[2])
        if curve_type == "POLY" or curve_type == "NURBS":
            point.co = (new_point_co[0], new_point_co[1], new_point_co[2], point.co[3])
    return curve_obj

def translate_curve_point(curve_obj, idx, trans_vector):
    curve_type = curve_obj.data.splines[0].type
    points = []
    if curve_type == "BEZIER":
        points = curve_obj.data.splines[0].bezier_points
    if curve_type == "POLY" or curve_type == "NURBS":
        points = curve_obj.data.splines[0].points
    for i in range(len(points)):
        if i == idx:
            points[i].co += trans_vector
            break

# Bezier curve!
# https://blender.stackexchange.com/questions/688/getting-the-list-of-points-that-describe-a-curve-without-converting-to-mesh
# https://gist.github.com/zeffii/5724956
def instance_on_bezier_curve(curve_obj, base_obj=None, scale_range=[1,2], frame_start=0, frame_end=1, collection_name=None):
    bezier_points = curve_obj.data.splines[0].bezier_points
    vert_per_segment = 15 # curve_obj.data.splines[0].resolution_u + 1
    leaf_scale_decrease = [0.5, 0.8]
    n_segments = len(bezier_points) - 1 # -1 if it is not cyclic!
    for i in range(n_segments):
        i_next = (i + 1) % len(bezier_points)
        point1 = bezier_points[i].co
        handle1 = bezier_points[i].handle_right
        handle2 = bezier_points[i_next].handle_left
        point2 = bezier_points[i_next].co
        interp_points = mathutils.geometry.interpolate_bezier(point1, handle1, handle2, point2, vert_per_segment)
        curve_tangent = mathutils.Vector(point1 - point2).normalized()
        curve_b, curve_n = pixar_onb(curve_tangent)
        for interp_point in interp_points:
            if mathutils.noise.random() > 0.7:
                continue
            leaf_rot_perturb = mathutils.Quaternion(curve_tangent, mathutils.noise.random()*3.14*2)
            curve_n.rotate(leaf_rot_perturb)
            inst = create_instance(base_obj, pos_vec=interp_point, rot_vec=curve_n, collection_name=collection_name)
            # Leaf scale animation must start later if leaf is further from root.
            inst.scale = (0,0,0)
            inst.keyframe_insert(data_path="scale", frame=frame_start) # from frame_start to adaptive_frame_start keep scale 0 and then start to increase scale.
            adaptive_frame_start = lerp(i / n_segments, frame_start, frame_end)
            inst.keyframe_insert(data_path="scale", frame=adaptive_frame_start)
            # Leaf scale must decrease with distance from root.
            scale_decrease = lerp(1 - i / n_segments, leaf_scale_decrease[0], leaf_scale_decrease[1])
            leaf_scale = lerp(mathutils.noise.random(), scale_range[0], scale_range[1]) * scale_decrease
            leaf_scale_vector = mathutils.Vector((leaf_scale, leaf_scale, leaf_scale))
            inst.scale = leaf_scale_vector
            inst.keyframe_insert(data_path="scale", frame=frame_end)
            #set_animation_fcurve(inst, "LINEAR")

#
# Constraint utils.
#

def instance_on_path(base_obj, scale_range, curve_path, n_inst, frame_start, frame_end, collection_name=None):
    dt = 1.0 / n_inst
    t = 0
    for i in range(n_inst):
        inst = create_instance(base_obj, rot_vec=mathutils.Vector((3,3,0)), collection_name=collection_name)
        inst.scale = (0,0,0)
        adaptive_frame_start = lerp(t, frame_start, frame_end)
        inst.keyframe_insert(data_path="scale", frame=adaptive_frame_start)
        leaf_scale = lerp(mathutils.noise.random(), scale_range[0], scale_range[1])
        leaf_scale_vector = mathutils.Vector((leaf_scale, leaf_scale, leaf_scale))
        inst.scale = leaf_scale_vector
        inst.keyframe_insert(data_path="scale", frame=frame_end)
        follow_path_constraint = inst.constraints.new("FOLLOW_PATH")
        follow_path_constraint.target = curve_path
        follow_path_constraint.use_fixed_location = True
        follow_path_constraint.offset_factor = t
        #set_animation_fcurve(inst, "LINEAR")
        t += dt

# Add and animate follow path constraint.
def follow_path(base_obj, path, offset_start, offset_end, frame_start, frame_end):
    follow_path_constraint = base_obj.constraints.new("FOLLOW_PATH")
    follow_path_constraint.target = path
    follow_path_constraint.use_fixed_location = True
    follow_path_constraint.offset_factor = offset_start
    follow_path_constraint.keyframe_insert(data_path="offset_factor", frame=frame_start)
    follow_path_constraint.offset_factor = offset_end
    follow_path_constraint.keyframe_insert(data_path="offset_factor", frame=frame_end)

#
# Animation utils.
#

# https://behreajj.medium.com/scripting-curves-in-blender-with-python-c487097efd13
def set_animation_fcurve(base_object, option='ELASTIC'):
    fcurves = base_object.data.animation_data.action.fcurves
    for fcurve in fcurves:
        for kf in fcurve.keyframe_points:
            # Options: ['CONSTANT', 'LINEAR', 'BEZIER', 'SINE',
            # 'QUAD', 'CUBIC', 'QUART', 'QUINT', 'EXPO', 'CIRC',
            # 'BACK', 'BOUNCE', 'ELASTIC']
            kf.interpolation = option
            # Options: ['AUTO', 'EASE_IN', 'EASE_OUT', 'EASE_IN_OUT']
            kf.easing = 'AUTO'

def animate_curve_bevel_growth(curve, frame_start, frame_end, growth_factor_end):
    curve.data.bevel_factor_end = 0
    curve.data.bevel_factor_start = 0
    curve.data.keyframe_insert(data_path="bevel_factor_end", frame=frame_start)
    curve.data.keyframe_insert(data_path="bevel_factor_start", frame=frame_start)
    curve.data.bevel_factor_start = growth_factor_end
    curve.data.keyframe_insert(data_path="bevel_factor_start", frame=frame_end)

def animate_curve_bevel_thickness(curve, frame_start, frame_end, bevel_min, bevel_max):
    curve.data.bevel_depth = 0.0
    curve.data.keyframe_insert(data_path="bevel_depth", frame=frame_start)
    curve.data.bevel_depth = lerp(mathutils.noise.random(), bevel_min, bevel_max)
    curve.data.keyframe_insert(data_path="bevel_depth", frame=frame_end)


def main():

    # Get sketch.
    sketch_curves = []
    for obj in bpy.data.collections["sketch_curves"].all_objects:
        sketch_curves.append(obj)

    # Get taper object.
    taper_obj = bpy.data.collections["taper_object"].all_objects[0]
    # Taper object is NURBS path and it has 5 points determining the shape of curve.
    taper_idx0_y_range = [-1, -0.5]
    taper_idx1_y_range = [-0.5, -0.25]
    taper_idx2_y_range = [-0.25, -0.12]
    taper_idx3_y_range = [-0.12, -0.05]
    taper_idx4_y_range = [-0.12, -0.05]

    # Get attractor.
    attractor = bpy.data.collections["attractor"].all_objects[0]
    attractor_scale_range = [0.1, 0.2]

    # Get firefly.
    firefly = bpy.data.collections["firefly"].all_objects[0]
    firefly_scale_range = [0.03, 0.05]

    # Get leaf.
    leaf = bpy.data.collections["leaf"].all_objects[0]
    leaf_scale_range = [0.15, 0.25]

    # Get fruit cluster.
    fruit_cluster = bpy.data.collections["fruit_cluster"].all_objects[0]
    fruit_cluster_scale_range = [0.05, 0.1]

    # Create perturbed copies of sketch and animate growth.
    n_copies_per_sketch = 5
    sketch_curve_copy_frame_start = 0
    sketch_curve_copy_frame_end = 400
    main_spline_bevel_range = [0.01, 0.02]
    sketch_copy_perturb_scale_range = [1.0, 1.2]
    sketch_copy_perturb_strength_range = [1.0, 1.1]
    sketch_copy_n_octaves = 2
    sketch_copy_perturb_amplitude_range = [0.5, 1.0]
    sketch_copy_perturb_frequency_range = [0.5, 1.0]
    sketch_copy_bevel_range = [0.01, 0.03]
    sketch_copy_growth_range = [0.7, 1.0]
    for sketch_curve in sketch_curves:
        # Sketch will be the largest spline.
        #animate_curve_bevel_growth(sketch_curve, sketch_curve_copy_frame_start, sketch_curve_copy_frame_end, 1.0)
        #animate_curve_bevel_thickness(sketch_curve, sketch_curve_copy_frame_start, sketch_curve_copy_frame_end, main_spline_bevel_range[0], main_spline_bevel_range[1])
        # Create attractor for sketch curve and add/animate follow path constraint.
        attractor_scale = lerp(mathutils.noise.random(), attractor_scale_range[0], attractor_scale_range[1])
        attractor_copy = create_instance(attractor, collection_name="attractor_copies")
        attractor_copy.scale = (attractor_scale, attractor_scale, attractor_scale)
        follow_path(attractor_copy, sketch_curve, 0.0, 1.0, sketch_curve_copy_frame_start, sketch_curve_copy_frame_end)
        #set_animation_fcurve(attractor_copy, "LINEAR")
        # Create sketch curve copies.
        for i in range(n_copies_per_sketch):
            # Create copy from sketch curve.
            sketch_copy = copy_object(sketch_curve, "sketch_curves_copies")
            # Perturb points of copied spline.
            sketch_copy_perturbed = perturb_curve(
                curve_obj=sketch_copy, 
                perturb_scale=lerp(mathutils.noise.random(), sketch_copy_perturb_scale_range[0], sketch_copy_perturb_scale_range[1]), 
                perturb_strength=lerp(mathutils.noise.random(), sketch_copy_perturb_strength_range[0], sketch_copy_perturb_strength_range[1]), 
                n_octaves=sketch_copy_n_octaves, 
                amplitude_scale=lerp(mathutils.noise.random(), sketch_copy_perturb_amplitude_range[0], sketch_copy_perturb_amplitude_range[1]), 
                frequency_scale=lerp(mathutils.noise.random(), sketch_copy_perturb_frequency_range[0], sketch_copy_perturb_frequency_range[1]))
            # Animate copy sketch spline growth.
            growth_factor_end = lerp(mathutils.noise.random(), sketch_copy_growth_range[0], sketch_copy_growth_range[1])
            animate_curve_bevel_growth(sketch_copy_perturbed, sketch_curve_copy_frame_start, sketch_curve_copy_frame_end, growth_factor_end)
            set_animation_fcurve(sketch_copy_perturbed, "LINEAR")
            # Add leafs to longer branches.
            if growth_factor_end > 0.0:
                instance_on_bezier_curve(curve_obj=sketch_copy_perturbed, base_obj=leaf, scale_range=leaf_scale_range, frame_start=sketch_curve_copy_frame_start, frame_end=sketch_curve_copy_frame_end, collection_name="leaf_copies")
            # Add fruits.
            if growth_factor_end > 0.9:
                instance_on_path(base_obj=fruit_cluster, scale_range=fruit_cluster_scale_range, curve_path=sketch_copy_perturbed, n_inst=7, frame_start=sketch_curve_copy_frame_start, frame_end=sketch_curve_copy_frame_end, collection_name="fruit_cluster_copies")
            # Add firefly at the top of copy spline.
            firefly_copy = create_instance(firefly, collection_name="firefly_copy")
            firefly_scale = lerp(mathutils.noise.random(), firefly_scale_range[0], firefly_scale_range[1])
            firefly_copy.scale = (firefly_scale, firefly_scale, firefly_scale)
            follow_path(firefly_copy, sketch_copy_perturbed, 0.0, growth_factor_end, sketch_curve_copy_frame_start, sketch_curve_copy_frame_end)
            #set_animation_fcurve(firefly_copy, "LINEAR")
            # Add bevel and animate.
            animate_curve_bevel_thickness(sketch_copy_perturbed, sketch_curve_copy_frame_start, sketch_curve_copy_frame_end, sketch_copy_bevel_range[0], sketch_copy_bevel_range[1])
            # Add taper object and configure its ponts.
            taper_copy = copy_object(taper_obj, "taper_object_copies")
            translate_curve_point(curve_obj=taper_copy, idx=0, trans_vector=mathutils.Vector((0.0, lerp(mathutils.noise.random(), taper_idx0_y_range[0], taper_idx0_y_range[1]), 0.0, 0.0)))
            translate_curve_point(curve_obj=taper_copy, idx=1, trans_vector=mathutils.Vector((0.0, lerp(mathutils.noise.random(), taper_idx1_y_range[0], taper_idx1_y_range[1]), 0.0, 0.0)))
            translate_curve_point(curve_obj=taper_copy, idx=2, trans_vector=mathutils.Vector((0.0, lerp(mathutils.noise.random(), taper_idx2_y_range[0], taper_idx2_y_range[1]), 0.0, 0.0)))
            translate_curve_point(curve_obj=taper_copy, idx=3, trans_vector=mathutils.Vector((0.0, lerp(mathutils.noise.random(), taper_idx3_y_range[0], taper_idx3_y_range[1]), 0.0, 0.0)))
            translate_curve_point(curve_obj=taper_copy, idx=4, trans_vector=mathutils.Vector((0.0, lerp(mathutils.noise.random(), taper_idx4_y_range[0], taper_idx4_y_range[1]), 0.0, 0.0)))
            #set_animation_fcurve(taper_copy, "LINEAR")
            sketch_copy_perturbed.data.taper_object = taper_copy


if __name__ == "__main__":
    main()