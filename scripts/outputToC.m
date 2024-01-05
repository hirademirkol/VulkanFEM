function outputToC(voxelized, fixed, loading, Ke0, name)
    global outPath_;
    out = permute(voxelized,[3 1 2]);
    writematrix(size(voxelized), strcat(outPath_, name, ".dat"), "Delimiter"," ");
    writematrix(out, strcat(outPath_, name, ".dat"),"WriteMode","append");
    
    writematrix(size(fixed, 1), strcat(outPath_, name, "_fixed.dat"));
    writematrix(fixed, strcat(outPath_, name, "_fixed.dat"),"WriteMode","append");

    writematrix(size(loading, 1), strcat(outPath_, name, "_loading.dat"));
    writematrix(loading, strcat(outPath_, name, "_loading.dat"), "Delimiter","\t","WriteMode","append");
    
	reOrdering = [10 11 12 1 2 3 7 8 9 4 5 6 22 23 24 13 14 15 19 20 21 16 17 18];
	Ke0 = Ke0(reOrdering, reOrdering);
    writematrix(Ke0, strcat(outPath_, name, "_Ke0.dat"));
end