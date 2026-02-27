function kep = kepVec2Struct(kepVec)
    kep.a = kepVec(1);
    kep.e = kepVec(2);
    kep.i = kepVec(3);
    kep.raan = kepVec(4);
    kep.w = kepVec(5);
    kep.theta = kepVec(6);
end