root -l <<EOF
gSystem->Load("libRooUnfold");
.L createResponseMatrixAngularity.C+
createResponseMatrixAngularity()
EOF
