import TDA
import Clustering
using Test

@testset "Clustering" begin

    # test knee point selection
    x = [0.470711, 0.215764, 0.17828, 0.085663, 0.0557541, 0.0279894, 0.0288017, 0.0178112, 0.0176353]
    @test last(TDA.knee(x)) == 4

    # test elbow selection method for clustering
    clusterfn = Clustering.kmeans
    clusterselectionfn = TDA.elbow
    Y = [-0.0790049 -0.10533 -0.16009 -0.17787 -0.199888 -0.235282 -0.244151 -0.283751 -0.25926 -0.214562 -0.2007 -0.176357 -0.115027 -0.0928281 -0.0955865;
          0.503826 0.510191 0.48794 0.488375 0.457246 0.47511 0.437912 0.436693 -0.411324 -0.418725 -0.405257 -0.439537 -0.463596 -0.438397 -0.457734]
    @test all(sort!(clusterselectionfn(clusterfn, Y)) .== vcat(fill(1,8), fill(2,7)))
    @test all(sort!(clusterselectionfn(clusterfn, Y, maxk=4)) .== vcat(fill(1,8), fill(2,7)))

end
