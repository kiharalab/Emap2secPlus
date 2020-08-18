#!/usr/bin/perl
if(@ARGV!=1){
 print "$0 [file]\n";
 exit;
}

($file)=@ARGV;

onetothree();

@B=readpdb("$file");

my $init=0;
my $chain_id="A";
my $residue_id="";
for($i=0;$i<@B;$i++){

  if($init==0&$B[$i][1] eq "CA "){
     $chain_id=$B[$i][7];
     printf(">$chain_id\n");
     $init=1;
  }
  if($chain_id ne $B[$i][7]& $B[$i][1] eq "CA "){
     $chain_id=$B[$i][7];
     printf("\n>$chain_id\n");
  }
 if($B[$i][1] eq "CA "&$B[$i][5] ne $residue_id){
  #print "$B[$i][0]= $B[$i][1] \n";
  printf("%s",$amin321{$B[$i][0]});
  $residue_id=$B[$i][5];
 }
}
print "\n";


#@key = sort { $hash{$a} <=> $hash{$b} || $a <=> $b} keys %hash;


sub firstfile{
my $cnt=0;
my @A;
open(IN,$_[0]) or die;
while(<IN>){
  next if(/^#/);
 chomp;
 my $item;
 @{$item}=split(/[\s\t]+/,$_);
 push @A, $item
}
close(IN);
return @A;
}
sub onetothree{
 %amin123 = ("W"=>"TRP","F"=>"PHE","Y"=>"TYR","L"=>"LEU","I"=>"ILE","V"=>"VAL","M"=>"MET","A"=>"ALA","G"=>"GLY","P"=>"PRO","C"=>"CYS","T"=>"THR","S"=>"SER","Q"=>"GLN","N"=>"ASN","E"=>"GLU","D"=>"ASP","H"=>"HIS","K"=>"LYS","R"=>"ARG");
 %amin321 = ("MSE"=>"M","HSD"=>"H","TRP"=>"W","PHE"=>"F","TYR"=>"Y","LEU"=>"L","ILE"=>"I","VAL"=>"V","MET"=>"M","ALA"=>"A","GLY"=>"G","PRO"=>"P","CYS"=>"C","THR"=>"T","SER"=>"S","GLN"=>"Q","ASN"=>"N","GLU"=>"E","ASP"=>"D","HIS"=>"H","LYS"=>"K","ARG"=>"R");
}
sub firstfile_line{
my $cnt=0;
my @A;
open(IN,$_[0]) or die;
while(<IN>){
  next if(/^#/);
 chomp;
 push @A, $_;
}
close(IN);
return @A;
}

sub readpdb{
 my $cnt=0;
 my @A;
 my ($x,$y,$z);
 my ($file)=@_;

 open(IN,$file) or die;
 while(<IN>){
  next unless(/^ATOM/);
  chomp;

  $x=substr($_,30,8);
  $y=substr($_,38,8);
  $z=substr($_,46,8);

  my $atm=substr($_,13,3);
  my $res=substr($_,17,3);
  my $rnum=substr($_,22,4);
  my $chain_id=substr($_,20,2);
  #my $m_tag=substr($_,17,9);
  my $m_tag=substr($_,13,13);

  my $item;
  @{$item}=($res,$atm,$x,$y,$z,$rnum,$m_tag,$chain_id);
  push @A, $item;
 }
 close(IN);
 return @A;
}


