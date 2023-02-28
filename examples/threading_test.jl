
function threaded(niters, s)
    a = zeros(niters)
    Threads.@threads for i = 1:niters
        A = rand(s, s)
        B = rand(s, s)
        a[i] = sum(A * B)
    end
end

function unthreaded(niters, s)
    a = zeros(niters)
    for i = 1:niters
        A = rand(s, s)
        B = rand(s, s)
        a[i] = sum(A * B)
    end
end

n = 10
s = 2000
@time threaded(n, s)
@time unthreaded(n, s)