function [v_rot] = rotation_vector_Rodrigues(vector, direction, delta)

v = vector;
u = direction;

v_rot = v*cos(delta) + cross(u, v)*sin(delta) + u.*dot(u, v).*(1-cos(delta));

end