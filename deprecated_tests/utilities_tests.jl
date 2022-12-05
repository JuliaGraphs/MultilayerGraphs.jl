struct test_struct
    a
    b
end

@test MultilayerGraphs.iscompletelyinitialized(test_struct(1, 2))
