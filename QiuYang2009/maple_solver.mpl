fname:=cat("degree",args[1],"solved.txt");
res:=solve(args[2]=0,a);
print(res);
save res, fname;
