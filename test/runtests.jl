using Test
import LinearAlgebra: dot
using RationalSimplex

@testset "RationalSimplex" begin

@testset "Test standard form" begin
      # Test case 1
      # min -10x - 12y - 12z
      # st     x +  2y +  2z <= 20
      #       2x +   y +  2z <= 20
      #       2x +  2y +   z <= 20
      #       x, y, z >= 0
      @testset "Test case 1 (optimal)" begin
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
      end

      # Test case 2
      # min -90x - y
      # st    2x + y <= 40
      #      20x + y <= 100
      #       2x     <= 3
      #        x,  y >= 0
      @testset "Test case 2 (optimal)" begin
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
      end

      # Test case 3
      # min x
      # st  x >= 3
      #     x <= 2
      @testset "Test case 3 (infeasible)" begin
            c = [1//1, 0//1, 0//1]
            A = [1//1 -1//1  0//1;
                 1//1  0//1  1//1]
            b = [3//1, 2//1]
            status, x = simplex(c, A, b)
            @test status == :Infeasible
      end

      # Test case 4
      # min -x
      # st  x >= 2
      @testset "Test Case 4 (unbounded)" begin
            c = [-1//1, 0//1]
            A = [ 1//1 -1//1;]
            b = [ 2//1]
            status, x = simplex(c, A, b)
            @test status == :Unbounded
      end

      # Check simplex detects nonneg rhs
      @testset "Check for nonnegative RHS" begin
            c = [-1//1, 0//1]
            A = [ 1//1 -1//1;]
            b = [-2//1]
            @test_throws AssertionError simplex(c, A, b)
      end
end

@testset "Test general form" begin
      # Test case 1
      # max  10x + 12y + 12z
      # st     x +  2y +  2z <= 20
      #       2x +   y +  2z <= 20
      #       2x +  2y +   z <= 20
      #       x, y, z >= 0
      @testset "Test case 1 (optimal)" begin
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
      end

      # Test case 2
      # max  90x + y
      # st    2x + y <= 40
      #      20x + y <= 100
      #       2x     <= 3
      #        x,  y >= 0
      @testset "Test case 2 (optimal)" begin
            c = [ 90//1,  1//1]
            A = [  2//1   1//1;
                  20//1   1//1;
                   2//1   0//1]
            b = [ 40//1, 100//1, 3//1]
            status, x = simplex(c, :Max, A, b, ['<','<','<'])
            @test status == :Optimal
            @test x[1] == 3//2
            @test x[2] == 37//1
      end

      # Test case 3
      # min x
      # st  x >= 3
      #     x <= 2
      @testset "Test case 3 (infeasible)" begin
            c = [1//1, 0//1]
            A = [1//1 0//1;
                 1//1 0//1]
            b = [3//1, 2//1]
            status, x = simplex(c, :Min, A, b, ['>','<'])
            @test status == :Infeasible
      end

      # Test case 4
      # min x
      # st   x >=  3
      #     -x >= -2
      @testset "Test case 4 (infeasible)" begin
            c = [1//1, 0//1]
            A = [ 1//1 0//1;
                 -1//1 0//1]
            b = [3//1, -2//1]
            status, x = simplex(c, :Min, A, b, ['>','>'])
            @test status == :Infeasible
      end
end

@testset "Kantorovich distances" begin
      # Problem description from 
      # http://stla.github.io/stlapblog/posts/KantorovichWithJulia.html
      c = [0//1, 1//1, 1//1, 1//1, 0//1, 1//1, 1//1, 1//1, 0//1]
      b = [1//7, 2//7, 4//7, 1//4, 1//4, 1//2]
      M = [1//1 1//1 1//1 0//1 0//1 0//1 0//1 0//1 0//1;
           0//1 0//1 0//1 1//1 1//1 1//1 0//1 0//1 0//1;
           0//1 0//1 0//1 0//1 0//1 0//1 1//1 1//1 1//1;
           1//1 0//1 0//1 1//1 0//1 0//1 1//1 0//1 0//1;
           0//1 1//1 0//1 0//1 1//1 0//1 0//1 1//1 0//1;
           0//1 0//1 1//1 0//1 0//1 1//1 0//1 0//1 1//1]

      status, x = simplex(c, :Min, M, b, ['=','=','=','=','=','='])
      @test dot(c,x) == 3//28
end

end  # RationalSimplex testset.