function vec = datetime2vec(t)
    vec = [year(t), month(t), day(t), hour(t), minute(t), second(t)];
end
