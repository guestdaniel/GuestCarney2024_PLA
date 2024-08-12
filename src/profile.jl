export test1 , test2

function test1()
    t = timevec(10.0, 10e3)
    τ = LogRange(1e-3, 1e0, 11)
    @time pea(t, τ)
end

function test2()
    t = timevec(10.0, 10e3)
    τ = LogRange(1e-3, 1e0, 11)
    @time pea2(t, τ)
end