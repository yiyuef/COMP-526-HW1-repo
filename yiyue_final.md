```julia
using Plots
using SymPy
using LinearAlgebra
using Printf
using Gridap
using GridapGmsh
using Gridap.Io
```


```julia
model123 = GmshDiscreteModel("/Users/apple/t1.msh")
# labels = get_face_labeling(model123)
```

    Info    : Reading '/Users/apple/t1.msh'...
    Info    : 11 entities
    Info    : 233 nodes
    Info    : 464 elements
    Info    : Done reading '/Users/apple/t1.msh'





    UnstructuredDiscreteModel()




```julia
# fn = "/Users/apple/model123.json"
# to_json_file(model123,fn)
```


```julia
# model123 = DiscreteModelFromFile("/Users/apple/model123.json")
# writevtk(model123,"/Users/apple/model123")
```


```julia
order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V0 = TestFESpace(model123,reffe;conformity=:H1,dirichlet_tags=["l1","l2","l3","l4","l5"])

```




    UnconstrainedFESpace()




```julia
boundary_cond = [500,400,300,200,100]
```




    5-element Vector{Int64}:
     500
     400
     300
     200
     100




```julia
# g(x) = 2.0
Ug = TrialFESpace(V0,boundary_cond)
```




    TrialFESpace()




```julia
degree = 2
Ω = Triangulation(model123)
dΩ = Measure(Ω,degree)
```




    GenericMeasure()




```julia
# neumanntags = ["l1","l2","l3","l4"]
# Γ = BoundaryTriangulation(model,tags=neumanntags)
# dΓ = Measure(Γ,degree)
```


```julia
# f(x) = 0


# a(u,v) = ∫( (∇(v)⋅∇(u)) - v⋅∇(u) )*dΩ
# b(v) =   ∫( 0*v )*dΩ

a(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ
b(v) = ∫( 0*v )*dΩ


```




    b (generic function with 1 method)




```julia
op = AffineFEOperator(a,b,Ug,V0)
uh = Gridap.Algebra.solve(op)
```




    SingleFieldFEFunction():
     num_cells: 414
     DomainStyle: ReferenceDomain()
     Triangulation: BodyFittedTriangulation()
     Triangulation id: 3200759391329012725




```julia
# ls = LUSolver()
# solver = LinearFESolver(ls)
```


```julia
# uh = Gridap.Algebra.solve(solver,op)
```


```julia
writevtk(Ω,"/Users/apple/results",cellfields=["uh"=>uh])
```




    (["/Users/apple/results.vtu"],)




```julia

```


```julia

```
