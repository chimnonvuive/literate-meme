function a = coriolis(r, omg, v_T)
  a = 2j*v_T*omg*exp(1j*angle(r));
end
