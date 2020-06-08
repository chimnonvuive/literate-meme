function v = vel(r, v_T, omg)
    v = v_T*exp(1j*angle(r)) + 1j*omg*r;
end