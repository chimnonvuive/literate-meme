function a = acc(r, v_T, omg, a_T, alp)
  a = a_T*exp(1j*angle(r)) + 1j*alp*r + coriolis(r,omg,v_T) + centripetal(r,omg);
end
