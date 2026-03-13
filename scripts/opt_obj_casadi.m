function [f,g]=opt_obj_casadi(f_function,x)
  [f,g] = f_function(x);
  f = full(f);
  g = full(g);
end