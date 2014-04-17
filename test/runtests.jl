using RationalSimplex
using Base.Test

# Test case 1
# min -10x - 12y - 12z
# st     x +  2y +  2z <= 20
#       2x +   y +  2z <= 20
#       2x +  2y +   z <= 20
#       x, y, z >= 0
c = [-10//1, -12//1, -12//1, 0//1, 0//1, 0//1]
A = [  1//1    2//1    2//1  1//1  0//1  0//1;
       2//1    1//1    2//1  0//1  1//1  0//1;
       2//1    2//1    1//1  0//1  0//1  1//1]
b = [20//1, 20//1, 20//1]
status, x = simplex(c, A, b)
@test status == :Optimal
@test x[1] == 4//1
@test x[2] == 4//1
@test x[3] == 4//1

# Test case 2
# min -90x - y
# st    2x + y <= 40
#      20x + y <= 100
#       2x     <= 3
#        x,  y >= 0
c = [-90//1, -1//1, 0//1, 0//1, 0//1]
A = [  2//1   1//1  1//1  0//1  0//1;
      20//1   1//1  0//1  1//1  0//1;
       2//1   0//1  0//1  0//1  1//1]
b = [ 40//1, 100//1, 3//1]
status, x = simplex(c, A, b)
@test status == :Optimal
@test x[1] == 3//2
@test x[2] == 37//1
@test x[3] == 0//1
@test x[4] == 33//1
@test x[5] == 0//1

# Test case 3
# min x
# st  x >= 3
#     x <= 2
c = [1//1, 0//1, 0//1]
A = [1//1 -1//1  0//1;
     1//1  0//1  1//1]
b = [3//1, 2//1]
status, x = simplex(c, A, b)
@test status == :Infeasible

# Test case 4
# min -x
# st  x >= 2
c = [-1//1, 0//1]
A = [ 1//1 -1//1;]
b = [ 2//1]
status, x = simplex(c, A, b)
@test status == :Unbounded

# Check simplex detects nonneg rhs
b = [-2//1]
@test_throws simplex(c, A, b)

##################################################

# Test case 1
# max  10x + 12y + 12z
# st     x +  2y +  2z <= 20
#       2x +   y +  2z <= 20
#       2x +  2y +   z <= 20
#       x, y, z >= 0
c = [ 10//1,  12//1,  12//1]
A = [  1//1    2//1    2//1;
       2//1    1//1    2//1;
       2//1    2//1    1//1]
b = [20//1, 20//1, 20//1]
status, x = simplex(c, :Max, A, b, ['<','<','<'])
@test status == :Optimal
@test x[1] == 4//1
@test x[2] == 4//1
@test x[3] == 4//1

# Test case 2
# max  90x + y
# st    2x + y <= 40
#      20x + y <= 100
#       2x     <= 3
#        x,  y >= 0
c = [ 90//1,  1//1]
A = [  2//1   1//1;
      20//1   1//1;
       2//1   0//1]
b = [ 40//1, 100//1, 3//1]
status, x = simplex(c, :Max, A, b, ['<','<','<'])
@test status == :Optimal
@test x[1] == 3//2
@test x[2] == 37//1

# Test case 3
# min x
# st  x >= 3
#     x <= 2
c = [1//1, 0//1]
A = [1//1 0//1;
     1//1 0//1]
b = [3//1, 2//1]
status, x = simplex(c, :Min, A, b, ['>','<'])
@test status == :Infeasible

# Test case 4
# min x
# st   x >=  3
#     -x >= -2
c = [1//1, 0//1]
A = [ 1//1 0//1;
     -1//1 0//1]
b = [3//1, -2//1]
status, x = simplex(c, :Min, A, b, ['>','>'])
@test status == :Infeasible