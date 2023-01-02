@test is_multiplex_interlayer(interlayer_multiplex_sg_mg)
@test is_multiplex_interlayer(interlayer_multiplex_sdg_vodg)

collect.(edges.(all_interlayers))
eltype.(all_interlayers)
edgetype.(all_interlayers)

for interlayer in all_interlayers
    @test all(has_edge.(Ref(interlayer), edges(interlayer)))
    @test all(has_vertex.(Ref(interlayer), mv_vertices(interlayer)))
    @test all(has_vertex.(Ref(interlayer), mv_vertices(interlayer)))
    mv = rand(collect(mv_vertices(interlayer)))
    mv_inneighbors(interlayer, mv)
    mv_outneighbors(interlayer, mv)
    adjacency_matrix(interlayer)
end

@debug ""

ne.(all_interlayers)
nv.(all_interlayers)

function _get_srcmv_dstmv_interlayer(interlayer::Interlayer)
    mvs = MultilayerGraphs.get_bare_mv.(collect(mv_vertices(interlayer)))

    src_mv_idx = findfirst(
        mv ->
            !isempty(
                setdiff(
                    Set(mvs),
                    Set(
                        vcat(
                            MultilayerGraphs.get_bare_mv.(mv_outneighbors(interlayer, mv)),
                            mv,
                            MultilayerGraphs.get_bare_mv.(mv_vertices(eval(mv.layer))),
                        ),
                    ),
                ),
            ),
        mvs,
    )

    src_mv = mvs[src_mv_idx]

    _collection = setdiff(
        Set(mvs),
        Set(
            vcat(
                MultilayerGraphs.get_bare_mv.(mv_outneighbors(interlayer, src_mv)),
                src_mv,
                MultilayerGraphs.get_bare_mv.(mv_vertices(eval(src_mv.layer))),
            ),
        ),
    )

    dst_mv = MultilayerGraphs.get_bare_mv(rand(_collection))

    return mvs, src_mv, dst_mv
end

interlayer = interlayer_sg_swg

mvs, src_mv, dst_mv = _get_srcmv_dstmv_interlayer(interlayer)

missing_edge = ME(src_mv, dst_mv)
# test uniform add_edge!
@test !has_edge(interlayer, missing_edge)
@test !has_edge(interlayer, MultilayerGraphs.reverse(missing_edge))
@test add_edge!(interlayer, missing_edge)
@test has_edge(interlayer, missing_edge)
@test has_edge(interlayer, MultilayerGraphs.reverse(missing_edge))
# test hybrid add_edge!
@test rem_edge!(interlayer, src_mv, dst_mv)
@test add_edge!(interlayer, src_mv, dst_mv)
@test has_edge(interlayer, src_mv, dst_mv)

interlayer = interlayer_swg_mg
mvs, src_mv, dst_mv = _get_srcmv_dstmv_interlayer(interlayer)
missing_edge = ME(src_mv, dst_mv, rand())
# test uniform add_edge!
@test !has_edge(interlayer, missing_edge)
@test add_edge!(interlayer, missing_edge)
@test has_edge(interlayer, missing_edge)
@test MultilayerGraphs.weights(interlayer)[
        interlayer.v_V_associations(src_mv), interlayer.v_V_associations(dst_mv)
    ] ==
    MultilayerGraphs.weights(interlayer)[
        interlayer.v_V_associations(dst_mv), interlayer.v_V_associations(src_mv)
    ] !=
    0.0
# test hybrid add_edge!
@test rem_edge!(interlayer, src_mv, dst_mv)
@test add_edge!(interlayer, src_mv, dst_mv, rand())
@test has_edge(interlayer, src_mv, dst_mv)
@test MultilayerGraphs.weights(interlayer)[
        interlayer.v_V_associations(src_mv), interlayer.v_V_associations(dst_mv)
    ] ==
    MultilayerGraphs.weights(interlayer)[
        interlayer.v_V_associations(dst_mv), interlayer.v_V_associations(src_mv)
    ] !=
    0.0

interlayer = interlayer_mg_vg
mvs, src_mv, dst_mv = _get_srcmv_dstmv_interlayer(interlayer)
missing_edge = ME(src_mv, dst_mv, (missing_metadata="hello",))
# test uniform add_edge!
@test !has_edge(interlayer, missing_edge)
@test add_edge!(interlayer, missing_edge)
@test has_edge(interlayer, missing_edge)

# test hybrid add_edge!
@test rem_edge!(interlayer, src_mv, dst_mv)
@test add_edge!(interlayer, src_mv, dst_mv, :missing_metadata, "hello")
@test has_edge(interlayer, src_mv, dst_mv)
@test get_prop(interlayer, src_mv, dst_mv, :missing_metadata) == "hello"

@test ne(interlayer_empty_sg_vg) == 0

interlayer = interlayer_sdg_swdg
mvs, src_mv, dst_mv = _get_srcmv_dstmv_interlayer(interlayer)
missing_edge = ME(src_mv, dst_mv)
# test uniform add_edge!
@test !has_edge(interlayer, missing_edge)
@test add_edge!(interlayer, missing_edge)
@test has_edge(interlayer, missing_edge)
# test hybrid add_edge!
@test rem_edge!(interlayer, src_mv, dst_mv)
@test add_edge!(interlayer, src_mv, dst_mv)
@test has_edge(interlayer, src_mv, dst_mv)

interlayer = interlayer_swdg_mdg
mvs, src_mv, dst_mv = _get_srcmv_dstmv_interlayer(interlayer)
missing_edge = ME(src_mv, dst_mv, rand())
# test uniform add_edge!
@test !has_edge(interlayer, missing_edge)
@test add_edge!(interlayer, missing_edge)
@test has_edge(interlayer, missing_edge)
@test MultilayerGraphs.weights(interlayer)[
    MultilayerGraphs.get_v(interlayer, src_mv), MultilayerGraphs.get_v(interlayer, dst_mv)
] != 0.0
# test hybrid add_edge!
@test rem_edge!(interlayer, src_mv, dst_mv)
@test add_edge!(interlayer, src_mv, dst_mv, rand())
@test has_edge(interlayer, src_mv, dst_mv)
@test MultilayerGraphs.weights(interlayer)[
    MultilayerGraphs.get_v(interlayer, src_mv), MultilayerGraphs.get_v(interlayer, dst_mv)
] != 0.0

interlayer = interlayer_mdg_vodg
mvs, src_mv, dst_mv = _get_srcmv_dstmv_interlayer(interlayer)
missing_edge = ME(src_mv, dst_mv, (weight=rand(),))
# test uniform add_edge!
@test !has_edge(interlayer, missing_edge)
@test add_edge!(interlayer, missing_edge)
@test has_edge(interlayer, missing_edge)
@test !isempty(get_metadata(interlayer, src_mv, dst_mv)) &&
    !isempty(get_metadata(interlayer, src_mv, dst_mv))
# test hybrid add_edge!
@test rem_edge!(interlayer, src_mv, dst_mv)
@test isempty(get_metadata(interlayer, src_mv, dst_mv)) &&
    isempty(get_metadata(interlayer, src_mv, dst_mv))
@test add_edge!(interlayer, src_mv, dst_mv, :weight, rand())
@test has_edge(interlayer, src_mv, dst_mv)
@test !isempty(get_metadata(interlayer, src_mv, dst_mv)) &&
    !isempty(get_metadata(interlayer, src_mv, dst_mv))

get_symmetric_interlayer.(all_interlayers)
