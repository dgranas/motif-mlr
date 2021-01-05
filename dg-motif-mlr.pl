#!/usr/bin/perl

#  Gurmukh Sahota
#  2008-07-10
#  Motif MLR with Missing Data
#  Beta 1.0

use FindBin qw($Bin);
use lib $FindBin::Bin;

my $USE_CONSTANT = 1;
my $GSL_MULTIFIT_C_PROGRAM = $Bin . '/GSL/multi_fit_r2';
my $INKSCAPE_PATH = $Bin . '/inkscape';
my $INKSCAPE_PROGRAM =  '/usr/bin/inkscape';

my $SAMPLER_ACTION_URL = "/cgi-bin/mkhan/sampling.pl";

my $USE_REGRESSION_VERSION = 2;
my $PRINT_REGRESSION_CSV = 0;
my $PRINT_COEFFICIENTS = 0;
my $PRINT_PVALUES      = 0;
my $ZOOM_LOGOS         = 0;
my $PRECISION = 2;
my $PALINDROMIC = 0;
my $REALLY_SMALL_NUMBER = 1e-200;
my $JOIN_CHAR = ',';
my $INDEX_OFFSET = 0;

my $LOGO_HEIGHT_CM = undef;
my $LOGO_WIDTH_CM  = undef;

use strict;
use Statistics::Regression;
use IPC::Open2;
use Data::Dumper;
use Memoize;
use Time::HiRes qw( usleep );
use SVGlogo qw( :ALL );
use MIME::Base64;

# memoize is a subroutine that allows perl to remember what the value of a function was
# given a set of input parameters (ie. memory of what input leads to what output).
memoize('return_coefficient_names', INSTALL => 'memoized_return_coefficient_names');
memoize('return_encoding_seq',      INSTALL => 'memoized_return_encoding_seq');


#
# Now we define a set of constants that are used later in the code
#
my $ALPHABET_SIZE = 4;
my $encoding_and_translation_h  = {
   "nucleotide" => {
                      "1" => {
                                "translation_h" => {
                                                      "seq2int" => {
                                                                      'A' => '0',
                                                                      'C' => '1',
                                                                      'G' => '2',
                                                                      'T' => '3'
                                                                   },
                                                      "int2seq" => {
                                                                      '0' => 'A',
                                                                      '1' => 'C',
                                                                      '2' => 'G',
                                                                      '3' => 'T'
                                                                   } 
                                                   },
                                "subencoding_hm"    => { 
                                                         "0" => [ 1, -1, -1],
                                                         "1" => [-1,  1, -1],
                                                         "2" => [-1, -1,  1],
                                                         "3" => [ 1,  1,  1]
                                                       },
                                "posname_a"            => ["W", "Y", "K"]
                             }
                   }
};


my @base_to_pseudobinary_translate_a;
my %base_to_int_h;
my %int_to_base_h;




#
# Below is part of the web/command-line interface capability (the output part)
# It checks to see if there are command line arguments and if so, it switches to command-line mode
# The "print" functionality is different between the two.
#
my $WEB = 1;
if (@ARGV) {
   $WEB = 0;
}

if ($WEB) {
   use CGI qw(:standard);
}
else {
   sub h2 {
       if ($WEB) {
         return CGI::h2($_[0]);
       }
       else {
         return '>>'   . $_[0] .   '<<' . "\n";
       }
   };
   sub h3 {
       if ($WEB) {
         return CGI::h3($_[0]);
       }
       else {
         return '>>>'  . $_[0] .  '<<<' . "\n";
       }
   };
   sub h4 {
       if ($WEB) {
         return CGI::h4($_[0]);
       }
       else {
         return '>>>>' . $_[0] . '<<<<' . "\n";
       }
   };
   sub b {
       if ($WEB) {
         return CGI::b($_[0]);
       }
       else {
         return $_[0];
       }
   };
   sub table {
       if ($WEB) {
         return CGI::table(@_);
       }
       else {
         return  (grep {! ref{$_} eq 'HASH'} @_);
       }
   };
   sub td {
       if ($WEB) {
         return CGI::td($_[0]);
       }
       else {
          return $_[0];
       }
   };
   sub Tr {
       if ($WEB) {
         return CGI::tr($_[0]);
       }
       else {
          return join("\t", @{$_[0]}). "\n";
       }
   };
}




# MAIN 
# This section of the code is simply i/o
# this is the second web/commandline interface
# in the end, both will call mlr_motif()

if ($WEB) {
#   if ($ENV{"HTTP_ACCEPT"} =~ /application\/xhtml\+xml/) {
#      print header('application/xhtml+xml');
#      print start_html(-title=>'Motif Model -- Multiple Linear Regression', -BGCOLOR=>'lightgrey'), "\n";
#   }
#   else {
      print header('text/html');
      if (param('Submit') eq 'Submit') {
         print start_html(-title=>'Motif Model -- Multiple Linear Regression', -BGCOLOR=>'lightgrey', -onLoad=>'moveTo(\'results\');'), "\n";
      }
      else {
         print start_html(-title=>'Motif Model -- Multiple Linear Regression', -BGCOLOR=>'lightgrey'), "\n";
      }
#   }
   print '<script type="text/javascript">
       var click_track = new Array();

       function moveTo(id)
       {
             if(document.getElementById && document.getElementById(id).scrollIntoView)
             {
                document.getElementById(id).scrollIntoView(true);
             };
             // will work in MSIE 5+ and Mozilla-based browsers. Not in Opera.
       }

       function setStyle(obj,style,value){
           getRef(obj).style[style]= value;
       }
       function getRef(obj){
           return (typeof obj == "string") ? document.getElementById(obj) : obj;
       }

      function selectOnClickTextfield(id, name, newValue)
      {
         if(click_track[id] == undefined) {
            text_box = document.getElementById(id);
            setStyle(text_box, \'color\', \'#000\');
            click_track[id] = true;
            text_box.value = "";
         }

         var radioObj = document.getElementsByName(name);
         if (!radioObj)
             return;
         var radioLength = radioObj.length;

         if (radioLength == undefined) {
           radioObj.checked = (radioObj.value == newValue.toString());
           return;
	 }
         for(var i = 0; i < radioLength; i++) {
	    radioObj[i].checked = false;
	    if(radioObj[i].value == newValue.toString()) {
		radioObj[i].checked = true;
 	    }
         }
      }
      function clearTextfieldOnClickOnce(id)
      {
         if(click_track[id] == undefined) {
            text_box = document.getElementById(id);
            setStyle(text_box, \'color\', \'#000\');
            click_track[id] = true;
            text_box.value = "";
         }
      }
      function setTextfieldOnClick(name, newValue)
      {
         var textObj = document.getElementsByName(name);
         if (!textObj)
             return;
         var textLength = textObj.length;

         if (textLength == undefined) {
           textObj.value = newValue.toString();
           return;
	 }
         for(var i = 0; i < textLength; i++) {
	    textObj[i].value = newValue.toString();
         }
      }
      </script>';

  my $default_data_message = "\n\nInput DNA Sequence followed by an Energy value (in kT or appropriate units)\nLower energies are considered to be tighter binders (Ex: Mnt binding protein)\n\n" . "GTGGACC 0\nGTGGCCC 0.41\nGTGGTCC 1\nGTGGACA 1.03\nGTGGACT 1.03\nGTGGCCA 1.44\nGTGGCCT 1.44\nGTGAACC 1.45\nGTGGACG 1.51\nGTGACCC 1.86\nGTGGCCG 1.92\nGTGGGCC 2\n";
  print
       '<center>', "\n",
       start_form( -name=>"motifmlr" ), "\n",
       "Enter sequence and data value separated by whitespace",  "\n", br,
       q(<a href="http://stormo.wustl.edu/cgi-bin/dgranas/motif_mlr.pl">Use this site if data includes methylated sequences</a>), # dave added this link
       p,  "\n",
       "";
       if ( (! defined param('mlr_data')) || (param('mlr_data') eq $default_data_message) ) { 
         print  textarea(-id => 'mlr_data', 
                         -name => 'mlr_data',
                         -default => "$default_data_message", 
                         -rows => 20, 
                         -cols => 80,
                         -style => 'color: #AAA;',
                         -onClick=>'clearTextfieldOnClickOnce(\'mlr_data\');'
                        );
       }
       else {
         print  textarea(-id => 'mlr_data', 
                         -name => 'mlr_data',
                         -rows => 20, 
                         -cols => 80,
                        );
       } 
       print 
       p, "\n",
       " Use file ", filefield('mlr_data_file', '', 50, 80), "\n",
       p, "\n",
       br, "\n",
       b("Fit to Models:" . "&nbsp;"x155),  "\n",
       br, "\n",
       b("Mono" . "&nbsp;"x7),
       radio_group( 
                    -name=>'mono',
                    -values=>['none', 'all', 'select'],
                    -default=>'all',
                    -labels=>{'none'=>'none ', 'all'=>'all ','select'=>''},
                    -onClick => 'setTextfieldOnClick(\'monoselect\', \'\');', 
                   ),
       textfield( 
                  -id=>'monoselect',
                  -name => 'monoselect', 
                  -default=>' 1   2   3   4   (Ex: index at 1, separate using space)',
                  -size => 60,
                  -maxlength => 80,
                  -style=> 'color: #AAA;',
                  -onClick => 'selectOnClickTextfield(\'monoselect\', \'mono\', \'select\');',
                ), 
       "\n", 
       br, "\n",
       b("Di" . "&nbsp;"x2),          
       radio_group( 
                    -name=>'di',
                    -values=>['none', 'all', 'adj', 'select'],
                    -default=>'none',
                    -labels=>{'none' => 'none ', 'all' => 'all', 'adj'=>'adj ', 'select'=>''},
                    -onClick => 'setTextfieldOnClick(\'diselect\', \'\');', 
                  ),
       textfield(
                  -id=>'diselect',
                  -name => 'diselect', 
                  -default=>'1,2   2,3   3,5  (Ex: index at 1, join with comma, separate using space)',
                  -size => 60,
                  -maxlength => 80,
                  -style=> 'color: #AAA;',
                  -onClick => 'selectOnClickTextfield(\'diselect\', \'di\', \'select\');',
                ), 
       br, "\n",
       b("Kmers" . "&nbsp;"x15),
       radio_group(
                    -name=>'kmer',
                    -values=>['none', 'select'],
                    -default=>'none',
                    -labels=>{'none' => 'none ' , 'select'=>''},
                    -onClick => 'setTextfieldOnClick(\'kmerselect\', \'\');', 
                   ),
       textfield(
                  -id => 'kmerselect',
                  -name => 'kmerselect', 
                  -default=>'1,2,6   1,2,3   (Ex: index at 1, join with comma, separate using space)',
                  -size => 60,
                  -maxlength => 80,
                  -style=> 'color: #AAA;',
                  -onClick => 'selectOnClickTextfield(\'kmerselect\', \'kmer\', \'select\');',
                ), 
       br, "\n",
       br, "\n",
       b("Symmetric/Palindromic Constraint". "&nbsp;"x5),
       radio_group(
                    -name=>'palindromic',
                    -values=>['off', 'on'],
                    -default=>'off',
                    -labels=>{'off' => 'off ' , 'on'=>'on '},
                   ),
       "&nbsp;"x5,
       b("Index Offset "),
       textfield(
                  -id => 'mlr_posoffset',
                  -name => 'mlr_posoffset', 
                  -default=>$INDEX_OFFSET,
                  -size => 3,
                  -maxlength => 3,
               ),
       "&nbsp;"x5,
       b("Display Decimal Precision "),
       textfield(
                  -id => 'precision',
                  -name => 'precision', 
                  -default=>$PRECISION,
                  -size => 3,
                  -maxlength => 3,
               ),
       "&nbsp;"x4,
       p, "\n",
       checkbox (
                  -name => 'coefficients',
                  -value => 'print',
                  -label => 'Print WYK Coefficients'
                ),
       "&nbsp;"x4,
       checkbox (
                  -name => 'pvalues',
                  -value => 'print',
                  -label => 'Print P-values'
                ),
       "&nbsp;"x4,
       checkbox (
                  -name => 'zoomlogos',
                  -value => 'yes',
                  -label => 'Zoom Logos',
                  -checked => 'true'
                ),
       "&nbsp;"x4,
       b("Logo Height (cm) "),
       textfield(
                  -id => 'logoheight',
                  -name => 'logoheight', 
                  -size => 4,
               ),
       "&nbsp;"x4,
       b("Logo Width (cm) "),
       textfield(
                  -id => 'logowidth',
                  -name => 'logowidth', 
                  -size => 4,
               ),
       "&nbsp;"x10,
       p, br, "\n",
       submit(-name=>'Submit', -value=>'Submit'), "\n",
       end_form, "\n",
       '</center>', "\n",
       hr, "\n";

   if (param('Submit') eq 'Submit') {
       print '<div id="results" />', "\n";

       my $data = param('mlr_data');
       my $data_filehandle = upload('mlr_data_file');
       if (defined $data_filehandle) {
         $data = '';
         while (<$data_filehandle>) {
            $data .= $_;
         }
       }
       #print "$data\n";
       my @models;
       foreach my $radio_param (qw(mono di kmer)) {
          my $radio_setting = param($radio_param);
          if ($radio_setting eq 'none') {
             # well, ignore it
          }
          elsif ($radio_setting eq 'select') {
             push @models, split /\s+/, param($radio_param . 'select');
          }
          else {
             push @models, ($radio_param . '_' . $radio_setting);
          }
       }

       my $coeff_print = param('coefficients');
       if ($coeff_print eq 'print') {
          $PRINT_COEFFICIENTS = 1;
       }
       my $pvalues_print = param('pvalues');
       if ($pvalues_print eq 'print') {
          $PRINT_COEFFICIENTS = 1;
          $PRINT_PVALUES = 1;
       }
       my $zoomlogos = param('zoomlogos');
       if ($zoomlogos eq 'yes') {
          $ZOOM_LOGOS = 1;
       }

       if (param('logoheight') ne '') {
          $LOGO_HEIGHT_CM = param('logoheight');
       }
       if (param('logowidth') ne '') {
          $LOGO_WIDTH_CM = param('logowidth');
       }

       if (param('mlr_posoffset') ne '') {
          $INDEX_OFFSET = param('mlr_posoffset');
       }

       my $palindromic_setting = param('palindromic');
       if ($palindromic_setting eq 'on') {
          $PALINDROMIC = 1;
       }

       $PRECISION = param('precision');


#       my @models = param('models');
       my $WM = mlr_motif($data, \@models);

       print hr;
   }
   print '</body></html>', "\n";
   print footer;
}
else {


   if ($ENV{"MLR_PRINT_COEFFICIENTS"} == 1) {
      $PRINT_COEFFICIENTS = 1;
   }
   if ($ENV{"MLR_PRINT_PVALUES"} == 1) {
      $PRINT_COEFFICIENTS = 1;
      $PRINT_PVALUES = 1;
   }
   if ($ENV{"MLR_PRECISION"}) {
      $PRECISION = $ENV{"MLR_PRECISION"};
   }

   my $datafile = shift @ARGV;
   my @models = @ARGV;
   open (IN, "<$datafile");
   my $data = '';
   while (<IN>) {
     $data .= $_;
   }
   close IN;
#   extend_mlr_encoding_info(3);
#   for my $nts (sort keys %{$$encoding_and_translation_h{"nucleotide"}}) {
#      print $nts, "...\n";
#      my $m = return_encoding_matrix($nts);
#      foreach my $row (@$m) {
#         print join(",", @$row), "\n";
#      }
#   }
#   exit(1);
   my $WM = mlr_motif($data, \@models);
}


# dg added subs

# Input: 
# - sequence
# - list of unmodified bases (ie ACGT)
# Output: True if sequence has modified base not 
sub seq_has_modified_base {
  my $seq = shift @_;
  my $dic_aref = shift @_;

  foreach my $char (split('', $seq)) {
        my $char_in_dic = 0; # set to 1 if char found in dict (ACGT)
        foreach my $dic_char (@$dic_aref) {
            #print "$char $dic_char\n";
            if ($char eq $dic_char) {
                $char_in_dic = 1;
                last; # break out of loop once found in dic (ie was ACGT)
            }
        }
        #if ($char_in_dic == 1) {
        #    print "was in dictionary\n";
        #}

        # if contained a non-ACGT char, return True that sequence has mod base
        if ($char_in_dic == 0) {
          return 1;
            # print "NOT in dictionary\n";
            # $has_mod_base = 1;
            # last;
        }
    }
    return 0; # no modified base was found
    # if ($has_mod_base == 1) {
    #     print "$s had a modified base";
    # }
}


#
# mlr_motif
# Input: sequence/energy data and models to run
# Output: motif model
# Purpose: This part of the code runs the error checking on the data input and then calls another subroutine
#          which will run the regression on each of the models selected
#
sub mlr_motif {

   my $data = shift @_;
   my $interactions_to_model = shift @_;

   my ($filtered_sequence_a, $filtered_response_a, $filtered_interactions_to_model_a, $err_value) = parse_input($data, $interactions_to_model);
	#print "@{$filtered_sequence_a}\n@{$filtered_response_a}\t@{$filtered_interactions_to_model_a}\n";
#   my $max_nonadjacency = 1;
#   foreach my $interaction_a (@$filtered_interactions_to_model_a) {
#      $max_nonadjacency = (scalar(@$interaction_a) > $max_nonadjacency) ? scalar(@$interaction_a) : $max_nonadjacency;
#   } 

#   &extend_mlr_encoding_info($max_nonadjacency);
   my $WM_model = mlr_regress($filtered_sequence_a, $filtered_response_a, $filtered_interactions_to_model_a);
   return $WM_model;
   
#  return (join(",", @sequence_array) . '-' . $response . '-' . join(",", @$interactions_to_model));
}


# 
# parse_input
# Input: sequence/energy data and models
# Output: sequence data, energy data, models and an error value
# Purpose: validate user input
#
sub parse_input {
   my $data_string = shift @_;
   my $models_a = shift @_;

   my $split_models_a;
   my @data_array = split /\n/, uc($data_string);
   
   my $seq_idx = 0;
   my $resp_idx = 1;
   my @sequence_array;
   my @response_array;
   foreach (@data_array) {
      $_ =~ s/^\s+//;
      $_ =~ s/\s+$//;
      my @elems = split /\s+/, $_;
      if ($elems[$seq_idx] !~ /^[ACTG]+$/) {
         print h3('Error: Rejected sequence: >' . $elems[$seq_idx]. '<' . ' ... Skipping line >' . $_ . '<');
         next;
      }
      if ($elems[$resp_idx] !~ /^\-?\d+\.?(\d+)*([eE]\-?\d+)?$/) {
         print h3('Error: Rejected response: >' . $elems[$resp_idx] . '<' . ' ... Skipping line>' . $_. '<' ); 
      }
      push @sequence_array, $elems[$seq_idx];
      push @response_array, $elems[$resp_idx];
   }
   
#   print h4('Info:: number of sequences (' . scalar(@filtered_sequence_a).')' ); 
#   print h4('Info:: number of sequences (' . scalar(@filtered_sequence_a).') which are (' . join(';', @filtered_sequence_a) . ')' ); 

   my $filtered_sequence_a;
   my $filtered_response_a;
   my $test_sequence_length = length($sequence_array[0]);
   for (my $idx=0; $idx < scalar(@sequence_array); $idx++) {
     
      my $len_seq = length($sequence_array[$idx]);
      if ($len_seq != $test_sequence_length) {
        print h3('Error: Sequence length for ' . $_ . '(' . $len_seq . ')'. ' does not equal ' . $test_sequence_length);
      }
      else {
         push @$filtered_sequence_a, $sequence_array[$idx];
         push @$filtered_response_a, $response_array[$idx];
      }
   }

   my $err_value = 0;
   if (scalar(@$filtered_sequence_a) == 0) {
     $err_value++;
   }

   if ($PALINDROMIC) {
      my $num_seq = scalar(@$filtered_sequence_a);
      for (my $idx=0; $idx < $num_seq; $idx++) {
         push @$filtered_sequence_a, reverse_complement($$filtered_sequence_a[$idx]);
         push @$filtered_response_a, $$filtered_response_a[$idx];
      }
   }
   foreach my $model (@$models_a) {
      if ($model eq 'mono_all') {
         for (my $i = 0; $i < $test_sequence_length; $i++) {
            push @$split_models_a, [$i+1];
         }
      }
      elsif ($model eq 'di_adj') {
         for (my $i = 0; $i < $test_sequence_length - 1; $i++) {
            push @$split_models_a, [$i+1, $i+2];
         }
      }
      elsif ($model eq 'di_all') {
         for (my $i = 0; $i < $test_sequence_length; $i++) {
            for (my $j = $i+1; $j < $test_sequence_length; $j++) {
               push @$split_models_a, [$i+1, $j+1];
            }
         }
      }
      else {
         if ($model !~  /^[0-9,]*[0-9]$/) {
            print h3('Error: Rejected model >' . $model . '< because of invalid specification');
         }
         else {
            push @$split_models_a, [map {$_ - $INDEX_OFFSET} (split /,/, $model)];
         }
      }
   }
   return $filtered_sequence_a, $filtered_response_a, $split_models_a, $err_value;
}

sub reverse_complement {
   my $string = shift @_;
   my $rc = reverse($string);
   $rc =~ tr/ACGT/TGCA/;
   return $rc;
}


sub return_coefficient_names {
   my @posmer = @_;
   #print "posmer = @posmer\n";
      
   my $len = scalar(@posmer);

   my $names;
   if ($len == 1) {
      foreach my $elem (@{$$encoding_and_translation_h{"nucleotide"}{"1"}{"posname_a"}}) {
         push @$names, $elem . ( $posmer[0] + $INDEX_OFFSET ) ;
      }
   }
   else {
      my @kmer = @posmer;
      my $last_pos = pop @kmer;

      my $kmer_code = memoized_return_coefficient_names(@kmer);
      my $last_code;
      foreach my $elem (@{$$encoding_and_translation_h{"nucleotide"}{"1"}{"posname_a"}}) {
         push @$last_code, $elem . ( $last_pos + $INDEX_OFFSET );
      }
    
      foreach my $elem (@$kmer_code) {
         foreach my $elem2 (@$last_code) {
            push @$names, ($elem . $elem2);
         }
      }

   }
   #print "names = @$names\n";
   return $names;
}

sub encode_subencoding_nmer {
   my $nmer = shift @_;
   #print "@{$nmer}";
   
   my $nts = scalar(@$nmer);

   my $nmer_sequence;
   for (my $i=0; $i<$nts; $i++) {
      $nmer_sequence .= $$encoding_and_translation_h{"nucleotide"}{"1"}{"translation_h"}{"int2seq"}{$$nmer[$i]};
   }
   #print "NMer_Seq: $nmer_sequence\n";
#   print STDERR "encode_subencoding_nmer :: $nts = $nmer_sequence\n";
   my $encoded; 
   if (exists $$encoding_and_translation_h{"nucleotide"}{$nts}{"translation_h"}{"seq2int"}{$nmer_sequence}) {
      $encoded = return_subencoding_nmer($nmer);
   }
   else {   
      my @kmer = @$nmer;
      my $last_nt = pop @kmer;

      my $kmer_code = encode_subencoding_nmer(\@kmer);
      my $nt_code = $$encoding_and_translation_h{"nucleotide"}{"1"}{"subencoding_hm"}{$last_nt};
    
      foreach my $elem (@$kmer_code) {
         foreach my $elem2 (@$nt_code) {
            push @$encoded, ($elem * $elem2);
         }
      }

      my $nts = scalar (@$nmer);
      my $row = base2dec($nmer, $ALPHABET_SIZE);
      $$encoding_and_translation_h{"nucleotide"}{$nts}{"subencoding_hm"}{$row} = [@$encoded];
      $$encoding_and_translation_h{"nucleotide"}{$nts}{"translation_h"}{"seq2int"}{$nmer_sequence} = $row;
      $$encoding_and_translation_h{"nucleotide"}{$nts}{"translation_h"}{"int2seq"}{$row} = $nmer_sequence;
   }

 # print "$$encoded[0] A $$encoded[1] A$$encoded[2] A $$encoded[3]\n";

 return $encoded;
}



sub print_hm_matrix
{
	my $hm_matrix = shift @_;

	
	foreach my $array_row (@$hm_matrix)
	{
		print "@$array_row", "\n";
		
	}

}







sub return_subencoding_nmer {
   my $nmer = shift @_;

   my $nts = scalar (@$nmer);
   my $row = base2dec(\@$nmer, $ALPHABET_SIZE);

   return $$encoding_and_translation_h{"nucleotide"}{$nts}{"subencoding_hm"}{$row};
}


sub return_encoding_seq {
   my $seq = shift @_;

   my $nts = length($seq);
   my $nmer;
   for (my $i=0; $i<$nts; $i++) {
      push @$nmer, $$encoding_and_translation_h{"nucleotide"}{"1"}{"translation_h"}{"seq2int"}{substr($seq, $i, 1)}
   }
   #print "@{$nmer}\n";
   my $encode_row;
#   print STDERR "TEST NMER:: @$nmer\n";

   if ($nts > 1) {
      my @rev_nmer = reverse @$nmer;
      my @powerset_kmer = sort {scalar(@$a) <=> scalar(@$b)} sub_powerset(@rev_nmer);
      foreach my $kmer (@powerset_kmer) {
         # due to hack in powerset code, we get a blank line instead of the full nmer
         next if (! scalar(@$kmer));
         my @rev_kmer = reverse (@$kmer);
#         print STDERR "TEST KMER:: @rev_kmer\n";
         push @$encode_row, @{encode_subencoding_nmer(\@rev_kmer)};
      }
   }
   push @$encode_row, @{encode_subencoding_nmer($nmer)};
   #print "Encode: @{$encode_row}\n";

   return $encode_row;
}

sub return_subencoding_hmatrix {
   my $nts = shift @_;

   return $$encoding_and_translation_h{"nucleotide"}{$nts}{"subencoding_hm"};
}

sub extend_seq2int {
   my $nts = shift @_;
   
   my @numeric_alphabet = sort {$a <=> $b} values %{$$encoding_and_translation_h{"nucleotide"}{"1"}{"translation_h"}{"seq2int"}};
   my @possible_nmers_a;
   for (my $i=1; $i<=$nts; $i++) {
      push @possible_nmers_a, [@numeric_alphabet];
   }
   my @nmers = recombine(@possible_nmers_a);
   foreach my $nmer (@nmers) {
      my $nmer_sequence = '';
      foreach my $elem (@$nmer) {
         $nmer_sequence .= $$encoding_and_translation_h{"nucleotide"}{"1"}{"translation_h"}{"int2seq"}{$elem};
      }
      next   if (exists $$encoding_and_translation_h{"nucleotide"}{$nts}{"translation_h"}{"seq2int"}{$nmer_sequence});
      my $idx = base2dec($nmer, $ALPHABET_SIZE);
      $$encoding_and_translation_h{"nucleotide"}{$nts}{"translation_h"}{"seq2int"}{$nmer_sequence} = $idx;
      $$encoding_and_translation_h{"nucleotide"}{$nts}{"translation_h"}{"int2seq"}{$idx} = $nmer_sequence;
   }
   return;
}

sub return_seq2int {
   my $nts = shift @_;
   
   return $$encoding_and_translation_h{"nucleotide"}{$nts}{"translation_h"}{"seq2int"};
}

sub return_int2seq {
   my $nts = shift @_;
   
   return $$encoding_and_translation_h{"nucleotide"}{$nts}{"translation_h"}{"int2seq"};
}


#
# mlr_regress()
# Input: array of sequences, array of responses and set of models
# Output: Hash of weight matrices (each hash key corresponds to a model)
# Purpose: Basically you can think of this as a large switch/case/if-then block that calls the individual model processing code
#          This and latter code is probably where attention needs to be focused for future development (there is alot of redundancy)
#
sub mlr_regress {
   my $sequences_a             = shift @_;
   my $response_a              = shift @_;
   my $interactions_to_model_a = shift @_;

   if (scalar(@$sequences_a) == 0) {
      print h3('No sequences in the file');
      return;
   }

   my $seq_len = length($$sequences_a[0]);

   my $WM_model = additive_nonadditive_mlr_model($sequences_a, $response_a, $interactions_to_model_a);
   #print "$WM_model\n";
#   print  Dumper($WM_model);
   if ($WEB) {
      print br,br, "\n";
      print start_form(-name=>'sample__all', -action => $SAMPLER_ACTION_URL), "\n";
      print submit (-name=>'Submit', -value=>'Sample using all PWMs');

      my $sample_data = '';
      foreach my $nts (sort keys %{$$WM_model{"PWM_full"}}) {
         my $pos2name_h = $$WM_model{"pos2name"}{$nts};
         my $pwm = $$WM_model{"PWM_full"}{$nts};
         $sample_data .= return_transpose_matrix($pwm, $pos2name_h);
      }
      print hidden(-name=> 'sample_data', -default=>[$sample_data]);
      print hidden(-name=> 'offset', -default=>[$$WM_model{"Intercept"}])   if ($USE_CONSTANT);
      print end_form;
   }

   my $info = [ 
                 ["R-squared", round($$WM_model{"R2"}, $PRECISION)], 
                 ["Intercept", round($$WM_model{"Intercept"}, $PRECISION)]
              ];
   if ($WEB) {
      print '<table>', "\n";
      print '<tr align="center">', "\n";
      print '<td>', "\n";
      web_print_table($info);
      print '</td>', "\n";
      print '<td>', "\n";
      web_print_OvE_plot($$WM_model{"Observed"}, $$WM_model{"Predicted"});
      print '</td>', "\n";
      print '</tr>', "\n";
      print '</table>', "\n";
      print br();
   }
   else {
      print_table($info);
   }
   if ($PRINT_COEFFICIENTS) {
      $info = [
                 ["Names",        (map{round($_, $PRECISION)} @{$$WM_model{"CoeffNames"}})],
                 ["Coefficients", (map{round($_, $PRECISION)} @{$$WM_model{"Coeffs"}})],
              ];
      if ($PRINT_PVALUES) {
         push @$info, ["P-values",     (map{round($_, $PRECISION)} @{$$WM_model{"Pvals"}})];
      }

      if ($WEB) {
         web_print_table($info);
      }
      else {
#         print "Coefficients:", "\t", join(",", @{$$WM_model{"Coeffs"}}), "\n";
#         if ($PRINT_PVALUES) {
#            print "P-values:",     "\t", join(",", @{$$WM_model{"Pvals"}}), "\n";
#         }
         print_table($info);
      }
   }
   my $iter = 0;
   foreach my $nts (sort keys %{$$WM_model{"PWM"}}) {
      $iter++;
#      print ">> ", $pwm_type, " :: ", $nts, "\n";
      my $pwm = $$WM_model{"PWM"}{$nts};

      my $pwm_type = 'PWM_' . $nts;
      my $pos2name_h = $$WM_model{"pos2name"}{$nts};

      if ($WEB) {
         print '<table>', "\n";
         print '<tr align="center">', "\n";
         print '<td>', "\n";
         my $form_name = "sample__" . $pwm_type;
         my $submit_value = 'Sample using ' . $pwm_type;
         print start_form(-name=>$form_name, -action => $SAMPLER_ACTION_URL), "\n";
         print submit (-name=>'Submit', -value=>$submit_value);

         my $sample_data .= return_transpose_matrix($pwm, $pos2name_h);
         print hidden(-name=> 'sample_data', -default=>[$sample_data]);
         print hidden(-name=> 'offset', -default=>[$$WM_model{"Intercept"}])   if ($USE_CONSTANT);
         print end_form;

         web_print_matrix($pwm, '<b>' . $pwm_type . '</b>', return_int2seq($nts), $pos2name_h);
         print '</td>', "\n";
         print '<td>', "\n";
         if ($iter == 1) {
            web_print_logo($pwm, $pos2name_h);
         }
         else {
            web_print_elogo($pwm, $pos2name_h);
         }
         print '</td>', "\n";
         print '</tr>', "\n";
         print '</table>', "\n";
         print br();
      }
      else {
         print_matrix($pwm, $pwm_type, return_int2seq($nts), $pos2name_h);
      }
   }
#      print h3('Test:: Completed Additive model ...');

   return;
}




# As you are reading this part of the code, I am going to define a set of functions 
# that will be used for the rest of the programs first and then go to each 
# subroutine (rather than in the order it is called).



# http://www.perlmonks.org/?node_id=197008
# hacked a little to simply return a set of combinations.
sub sub_powerset {
    my @list= @_;           # List of items to choose from
    my @pick= (0) x @list;  # Whether we want each item
    # $pick[$i] means include $list[$i] in results.
    # So @pick currently describes the empty subset.
    # Return a closure that, each time it is called, returns
    #   the next subset:

    my @final_a;
    while (1) {
        # Treat @pick as a base-2 number and increment it.
        # Note that @pick started as all 0s and we stop
        #   after it is all 1s so all cases get covered.
        # (See original node for handling the empty subset)

        # Start at least-significant bit, $pick[0]:
        my $i= 0;
        # Increment a bit. If the bit was already 1, then
        #   set it to 0 and continue to next bit:
        while( 1 < ++$pick[$i]  ) {
            $pick[$i]= 0;
            # If we've run out of bits, then we were at
            #   all 1s and so are done. Return empty list:
            return @final_a  if  $#pick < ++$i;
        }
        # The grep() below returns the indices for which
        #   $pick[$_] is not 0.  The @list[...] is an array
        #   slice that returns the list of elements of @list
        #   at the indices returned by grep.  That is, we
        #   return all items $list[$i] where $pick[$i] is
        #   not 0.  Same as:
        #     map { $pick[$_] ? $list[$_] : () } 0..$#list;
        push @final_a, [@list[ grep {!$pick[$_]} 0..$#pick ]];
    }
    return;
}


#---------------------------------------------------------------------------------------------------------------#
# sub recombine("array of array of character strings" @strings )                                                #
#                                                                                                               #
#       Input   : array of array of character elements.                                                         #
#       Output  : an array of an array of recombined elements. One from each of the outer array.                #
#       Purpose : To recombine each of the inner array elements with the elements of the others.                #
#                      One from each outer array. Makes atom filters in this script.                            #
#                                                                                                               #
#       How     : an array of current, an array of maxes.  You increment the last element until it maxes out.   #
#                      Then you run though the counter in reverse.  If the element is maxed.  Add one to the    #
#                      previous element.  If the first element is maxed then exit the sub because you are done  #
#---------------------------------------------------------------------------------------------------------------#

sub recombine
  {
  # the array of the strings
  my @strings = @_;

  # initialize the variables
  my $last_string = $#strings;     
  my @max_elements;
  my @current;
  my @recombination;
  my @element_of_recombination;
  my $count;
  my $element_count;
  
  # initalize the maximum number of elements for the inner arrays and sets the current to 0.
  foreach my $elements_of_string (@strings)
  {
      if ($#$elements_of_string == -1)
      {
	  push @max_elements, 0;
      }
      else
      {
	  push @max_elements, $#$elements_of_string;
      }
      push @current, 0;
  }
  # an infinite while loop which exits when the first inner array's counter maxes out.
  while (1)
    {
    # Run through each element of the last string  
    for ($element_count=0; $element_count<=$max_elements[$last_string];$element_count++)
      {
      # Go through each element of the outer array and build a list of the recombined elements
      my $temp_string_count = 0;
      foreach my $elements_of_string (@strings)
        {
        push @element_of_recombination, $$elements_of_string[$current[$temp_string_count]];  
        $temp_string_count++;
        }
      # push the list of recombined elements onto your return array of arrays.  
      push @recombination, [@element_of_recombination];
      @element_of_recombination = ();
      $current[$last_string]++;
      }
    # Do the "carries" of the current elements 
    for ($count=$last_string; $count>0; $count--)
      {
      # if the current count is greater than the max number of elements,
      #     add one to prev and set current to 0.  
      if ($current[$count] > $max_elements[$count])
        {
        $current[$count] = 0;
        $current[$count-1] += 1;
        }
      # Once you hit an element that does not need to be carried exit.
      else
        {last;}
      }
    # if your first element has maxed out then exit the script
    if ($current[0] > $max_elements[0])
      {return (@recombination);}
    }
  
  }



# 
# mlr()
# Input: names, y variables and X rows (ie. regression variables), as an array of y's and matrix of X respectively
# Output: an r2 fit and a theta (ie. regression coefficients)
# Purpose: to run the multiple linear regression
#          We can see that you can either use the GSL based code that Barrett Foat wrote, or the inbuilt perl module
#          Barrett's GSL C stuff is much faster and throws less errors than the PERL version, there is probably a better way to 
#          Incorporate it (like inline C), but I never got around to it.
#
sub mlr {
   my $n = shift @_;
   my $y = shift @_;
   my $X = shift @_;

   my $r2;
   my $theta;
   my $pval;
   my $predict;

   my $t_values;

   if ($PRINT_REGRESSION_CSV) {
      for (my $count=0; $count <@$y; $count++) {
         print join(",", ( $$y[$count], @{$$X[$count]} )), "\n";
      }
   }

   if ($USE_REGRESSION_VERSION == 0) {
      my $reg = Statistics::Regression->new( "sample regression", $n );
      for (my $count=0; $count <@$y; $count++) {
         $reg->include( $$y[$count], [@{$$X[$count]}] );
      }
      @$theta = $reg->theta();
      $r2 = $reg->rsq();
   }
   elsif ($USE_REGRESSION_VERSION == 1) {
      ($r2, $t_values, $theta) = barrett_gsl_multivariate_fit($y, $X);
   }
   elsif ($USE_REGRESSION_VERSION == 2) {
      ($r2, $theta, $pval, $predict) = use_R_mlr($y, $X);
   }

#   print "MLR => ", $r2, "\t", join(",", @$theta), "\n";
   return ($r2, $theta, $pval, $predict);
}


#
# Basically calls the C program that Barrett wrote and does the appropriate I/O for the mlr
#
sub barrett_gsl_multivariate_fit {
  my ($ra_y, $ra_X) = @_;

  #my $multi_fit = '/home/barrett/projects/LoopE/c/multi_fit_r2';
  my $multi_fit = $GSL_MULTIFIT_C_PROGRAM;

  #Recruit laborer
  local (*READ, *WRITE);
  my $pid = open2(\*READ, \*WRITE, $multi_fit);
  #open WRITE, ">test.pipe";
  my $n = @$ra_X;
  my $p = @{$$ra_X[0]};

  print WRITE "$n,$p\n";
  for (my $i = 0; $i < $n; $i++){
    print WRITE "$$ra_y[$i]\n";
    for (my $j = 0; $j < $p; $j++){
      print WRITE "$$ra_X[$i][$j]\n";
    }
  }

  my @Fs;
  my @ts;
#=pod
  my $r2 = <READ>;
  chomp $r2;
  while (<READ>){
    chomp;
    my ($F, $t) = split /\t/, $_;
    push @Fs, $F;
    push @ts, $t;
  }
  close WRITE;
  close READ;
  waitpid($pid, 0); #You're fired!
#=cut
#  close WRITE; die;

  return ($r2, \@ts, \@Fs);
}


sub use_R_mlr {
   my ($y, $X) = @_;

   my $R_EXECUTABLE = 'R --vanilla -q';
#   $SIG{PIPE} = 'IGNORE';
   #Recruit laborer
   local (*READ, *WRITE);
   my $pid = open2(\*READ, \*WRITE, $R_EXECUTABLE);
   #
   
   print WRITE 'wyk=read.table(stdin(),sep=",")', "\n";
   my $trash = <READ>;
   my $examples = scalar(@$y);
   my $row_string = '1:' . $examples;
   my $Xlen = scalar(@{$$X[0]});
   $Xlen-- if ($USE_CONSTANT);
   for (my $i=0; $i<$examples; $i++) {
      my @x = @{$$X[$i]};
#      print STDERR "$i / $examples ...\n";
#      print STDERR @x, "\n";
      if ($USE_CONSTANT) {
         shift @x;
      }
      print WRITE join(",", ($$y[$i], @x)) .  "\n";
      $trash = <READ>;
   }
   print WRITE "\n";
   $trash = <READ>;
   my $dependants = '';
   for (my $i=2; $i<=$Xlen+1; $i++) {
      $dependants .= 'wyk[' . $row_string . ',' . $i . ']+';
   }
   $dependants =~ s/\+$//;

   print WRITE 'attach(wyk)', "\n";
   $trash = <READ>;
   print WRITE 'mlr=lm(wyk[' . $row_string . ',1]' . '~' . $dependants .')', "\n";
   $trash = <READ>;
   print WRITE 'summary(mlr)$r.squared', "\n";
 
  my $r2 = 0;
   my $line = <READ>;
   if ($line =~ /^>\s+summary\(mlr\)\$r\.squared/) {
      my $temp = <READ>;
      chomp $temp;
      my @elems =split /\s+/, $temp;
      $r2 = $elems[1];
   }

#   print WRITE 'coeff=mlr$coef', "\n";
#   my $trash = <READ>;
   print WRITE 'write.table(summary(mlr)$coefficients)', "\n";
   print WRITE 'write.table(predict(mlr))', "\n";
   print WRITE 'q()', "\n";

   my $intercept = 0;
   my $ipval = 0;
   my $coeff;
   my $pval;
   my $predict;
   for (my $i=2; $i<=$Xlen+1; $i++) {
      $$coeff[$i-2] = '-';
      $$pval[$i-2]  = '-';
   }

   while (my $line = <READ>){
      chomp $line;
      if ($line =~ /^>\s+write.table\(summary\(mlr\)\$coefficients\)/) {
         <READ>;
         while (my $line2 = <READ>) {
            chomp $line2;

            if ($line2 =~ /^\>/) {
               $line = $line2;
               last;
            }
            if ($line2 =~ /^\"\(Intercept\)"\s+(.*)\s+(.*)\s+(.*)\s+(.*)$/) {
               $intercept = $1;
            }
            elsif ($line2 =~ /^\"wyk\[$row_string\, (\d+)\]"\s+(.*)\s+(.*)\s+(.*)\s+(.*)$/) {
               my $idx = $1-2;
               my $val = $2;
               $$coeff[$idx] = $val;
               $$pval[$idx] = $5;
            }
            else {
               print STDERR "ERROR :: R ... cannot parse $line2 in write.table coefficients parsing\n";
            }
         }
      }
      if ($line =~ /^>\s+write.table\(predict\(mlr\)\)/) {
         <READ>;
         while (my $line2 = <READ>) {
            chomp $line2;

            if ($line2 =~ /^\>/) {
               $line = $line2;
               last;
            }
            if ($line2 =~ /^\"(\d+)\"\s+(.*)$/) {
               my $idx = $1-1;
               my $val = $2;
               $$predict[$idx] = $val;
            }
            else {
               print STDERR "ERROR :: R ... cannot parse $line2 in write.table prediction parsing\n";
            }
         }
      }
   }

   #print WRITE 'q()', "\n";
   #$trash = <READ>;
   unshift @$coeff, $intercept    if ($USE_CONSTANT);

   close WRITE;
   close READ;

   waitpid($pid, 0); #You're fired!
 
   return ($r2, $coeff, $pval, $predict);
}

#
# The next set of functions are simply for converting between different base numbers (ie. base 4 -> base 10).
# The reason for this is that I wanted to go from the nucleotide coding (ie. base 3/4) to an array index (base 10).
#
sub logb {
   my $BASE = shift @_;
   my $value = shift @_;

   return (log($value) / log($BASE));
}

sub base2dec {
   my $arr = shift @_;
   my $BASE = shift @_;

   my $dec = 0;
   for (my $count=0; $count<scalar(@$arr); $count++) {
      $dec += $$arr[$count] * ($BASE ** ($#$arr-$count));
   }

   return $dec;
}

sub dec2base {
   my $dec = shift @_;
   my $BASE = shift @_;
   my $precision = shift @_;

   my @arr;
  
   my $num = $dec;
   for (my $count=($precision-1); $count>=0; $count--) {
      my $place = ($BASE ** $count);
      my $div = int($num / $place);
      push @arr, $div;
      $num -= ($div * $place);
   }
   return [@arr];
}

################################### round ##############################
#                                                                      #
# Input   : Variable, log(10^x) power (1 for 0.1 ...)                  #
# Output  : Rounded Variable to the 10^ power (only positive integers) #
# Purpose : Rounding variables                                         #
########################################################################

sub round 
   {
       my $rounding_var = shift @_;
       my $power = shift @_;
       
       # failsafes
       return $rounding_var    if ($rounding_var eq '-');
       return $rounding_var    if ($power < 0);
       return $rounding_var    if ($rounding_var !~ /^\-??[0-9]*\.??[0-9]+(e\-??\d+)??$/);

       my $round_offset = 0.5;
       if ($rounding_var < 0) {
          $round_offset = -0.5;
       }
       my $string = '%.' . $power . 'f';
       # This sprintf statement was the only way I could keep the 30.0 from becoming 30.  
       return sprintf($string, (((int(($rounding_var * (10 ** $power)) + $round_offset))) / (10 ** $power)));
   }



#
# The next 4 subroutines are basically the separate types of models that are currently available:
#    additive_missing_mlr_model
#    additiveANDadjacentdi_missing_mlr_model
#    additive_mlr_model
#    additiveANDadjacentdi_mlr_model
#    They all have the same basic structure (the most recent and "best" designed code is probably from the missing models):
#        a note of caution, the newer missing models probably would not be able to be used with the Perl statistics subroutines
#        as they never define the regression variable names (@regression_var_names).  This is probably an easy fix though.
#
#    The flow is as follows:
#          1. Get Input data
#          2. create regression variables (ie. X, y and names)
#                this requires the coding scheme and is where alot of the design for the "list-based" non-additivity will be
#          3. missing data removal/analysis
#                this code is pretty intelligent for what it does now, but will probably also have to be redesigned to handle the non-additivity
#          4. run MLR
#                you could update this part to only use the C code (and hopefully better integrate it with the PERL)
#          5. parse out the regression results into a PWM
#                this looks ugly, but works for now.  you will probably have to do alot of printing to figure out exactly how it works.
#          6. return



sub additive_nonadditive_mlr_model {
   my $sequences_a              = shift @_;
   my $response_a               = shift @_;
   my $interactions_to_model_a  = shift @_;

   my $WM;
   my @regression_var_names;
#   print h2('Test:: In additive_mlr_model');
   my $sequence_length = length ($$sequences_a[0]);

   # Create regression object
   my @global_y;
   for (my $count=0; $count<scalar (@$sequences_a); $count++) {
      push @global_y, $$response_a[$count];
   }

#   print STDERR "ENCODING ....\n";
   my $X_h;
   my $count_nts_interaction_h;
   my $interaction_nts_h;

   foreach my $interaction (@{$interactions_to_model_a}) {
        #print "@$interaction\n";
      	my $nts = scalar (@$interaction);
        #print "$nts\n";
        $$count_nts_interaction_h{$nts}++;
        push @{$$interaction_nts_h{$nts}}, $interaction;
        for (my $count=0; $count<scalar(@$sequences_a); $count++) {
        	my @encode_a;
         	my $seq_substr;
         	foreach my $seq_pos (@$interaction) {
                	#print "$seq_pos\n";
         		$seq_substr .= substr($$sequences_a[$count], ($seq_pos-1), 1);
         		#print "$seq_substr\n";
         	}
        	#print $seq_substr, ">>", @{return_encoding_seq($seq_substr)}, "\n";
         	push @{$$X_h{$nts}[$count]}, @{memoized_return_encoding_seq($seq_substr)};
        }
   }

   my @nts_a = sort {$a <=> $b} keys %$X_h;
  
   for (my $nts=2; $nts<=$nts_a[-1]; $nts++) {
      my $seq2int_h = return_seq2int($nts);
      my $num_mapped_seqs = scalar(keys %$seq2int_h);
      my $num_potential_seqs = $ALPHABET_SIZE ** $nts;
      if ($num_mapped_seqs < $num_potential_seqs) {
         extend_seq2int($nts);
      }
   }
   my $superX;
   if ($USE_CONSTANT) {
      for (my $count=0; $count<scalar (@$sequences_a); $count++) {
         push @{$$superX[$count]}, 1.0;
      }
   }

   foreach my $nts (@nts_a) {
     for (my $count=0; $count<scalar (@$sequences_a); $count++) {
         push @{$$superX[$count]}, @{$$X_h{$nts}[$count]};
      }
   }

   my $d=0;
   foreach my $elem (@{$$superX[0]}) {
      push @regression_var_names, ('d' . $d++);
   }

#   print STDERR "REGRESSING ....\n";
   #ADD USE CONSTANT TO NEW_X HERE.
   my ($r2, $thetas, $pval, $predict) = &mlr(\@regression_var_names, \@global_y, $superX);

   
#   print STDERR "PARSING ....\n";

#   if ($PRECISION > 0) {
#      my $test_zero = 10**(-1 * $PRECISION);
#      foreach (@$thetas) {
#         next       if ($_ eq '-');
#
#         $_ = 0     if ($_ < $test_zero);
#      }
#   }

   my $constant_intercept=0;
   if ($USE_CONSTANT) {
      $constant_intercept=shift @$thetas;
      $$WM{"Intercept"} = $constant_intercept;
   }
   $$WM{"R2"} = $r2;
   $$WM{"CoeffNames"};
   @{$$WM{"Coeffs"}}    = @$thetas;
   @{$$WM{"Pvals"}}     = @$pval;
   @{$$WM{"Predicted"}} = @$predict;
   @{$$WM{"Observed"}}  = @global_y;

#   print "CONSTANT_INTERCEPT:: ", $constant_intercept, "\n";

   my $seen_interactions_h;
   my $coeffName_h;
   my $WM_type = 'PWM';
   foreach my $base_nts (@nts_a) {

#      my $max_seq_len = $$count_nts_interaction_h{$nts};
      foreach my $int (@{$$interaction_nts_h{$base_nts}}) {
         my @rev_int = reverse @$int;
         my @powerset_posns = sort {scalar(@$a) <=> scalar(@$b)} sub_powerset(@rev_int);
         push @powerset_posns, [@rev_int];
         my @position_sets;
         foreach my $posns (@powerset_posns) {
            # due to hack in powerset code, we get a blank line instead of the full nmer
            next if (! scalar(@$posns));
            push @position_sets, [reverse(@$posns)];
         }
         foreach my $position_set (@position_sets) {
            my $nts = scalar (@$position_set);
            my $seen_seq_h;
            for (my $count=0; $count<scalar (@$sequences_a); $count++) {
               my $seq_substr;
               foreach my $seq_pos (@$position_set) {
                  $seq_substr .= substr($$sequences_a[$count], ($seq_pos-1), 1);
               }
               $$seen_seq_h{$seq_substr}++;
            }
            my $name = join($JOIN_CHAR, (map {$_ + $INDEX_OFFSET} @$position_set));
            $$seen_interactions_h{$name}++;

            my $seq_pos = -1;
            if (! (exists $$WM{"name2pos"}{$nts}{$name})) {
               $seq_pos = scalar (keys %{$$WM{"name2pos"}{$nts}});
               $$WM{"name2pos"}{$nts}{$name} = $seq_pos;
            }
            else {
               $seq_pos = $$WM{"name2pos"}{$nts}{$name};
            }

            my $subencoding_hm = return_subencoding_hmatrix($nts);
            my @subencoding_hm_keys = keys %$subencoding_hm;
            my $encode_line = scalar(@{$$subencoding_hm{$subencoding_hm_keys[0]}});
            my $seq2int_h = return_seq2int($nts);

            my $num_usable_cols = $encode_line;
            my @posset_coeff_names = @{memoized_return_coefficient_names(@$position_set)};
            my @use_thetas = splice(@$thetas,0,$num_usable_cols);
            for (my $idx=0; $idx<scalar(@use_thetas); $idx++) {
               if ($use_thetas[$idx] eq '-') {
                  $use_thetas[$idx] = 0;
               }
               $$coeffName_h{$posset_coeff_names[$idx]} += $use_thetas[$idx];
            }
            push @{$$WM{"CoeffNames"}}, @posset_coeff_names;
#         print STDERR "Testing:: Usable columns $num_usable_cols out of ", $encode_line,"\n";
            foreach my $seq (sort {$$seq2int_h{$a} <=> $$seq2int_h{$b}} keys %$seq2int_h) {
               my $row = $$seq2int_h{$seq};
               if (exists $$seen_seq_h{$seq}) {
                  my @use_encodings = @{$$subencoding_hm{$row}};
                  if (scalar(@use_encodings) != scalar(@use_thetas)) {
                     print STDERR "ERROR:: Number of dependant encoding variables (", scalar(@use_encodings), ") does not equal the number of thetas(", scalar(@use_thetas), ") that were pulled out\n";
                  }
                  else {
                     $$WM{$WM_type}{$nts}[$row][$seq_pos] += 0;
                     for (my $idx=0; $idx<scalar(@use_thetas); $idx++) {
                        $$WM{$WM_type}{$nts}[$row][$seq_pos] += $use_encodings[$idx]*$use_thetas[$idx];
                     }
                  }
               }
               else {
                  $$WM{$WM_type}{$nts}[$row][$seq_pos] = '-';
               }
            }
         }
      }
   }
   foreach my $nts (keys %{$$WM{"name2pos"}}) {
      foreach my $name (keys %{$$WM{"name2pos"}{$nts}}) {
         $$WM{"pos2name"}{$nts}{$$WM{"name2pos"}{$nts}{$name}} = $name;
      }
   }


   @$thetas = @{$$WM{"Coeffs"}};
   my $WM_type = 'PWM_full';
   foreach my $base_nts (@nts_a) {

#      my $max_seq_len = $$count_nts_interaction_h{$nts};
      foreach my $int (@{$$interaction_nts_h{$base_nts}}) {
         my $use_interaction = 1; 
         my $full_name = join($JOIN_CHAR, (map {$_ + $INDEX_OFFSET} @$int));
         my $full_nts = scalar(@$int);
         my $full_seq2int_h = return_seq2int($full_nts);
         $use_interaction = 0    if ($$seen_interactions_h{$full_name} > 1);
         my @rev_int = reverse @$int;
         my @powerset_posns = sort {scalar(@$a) <=> scalar(@$b)} sub_powerset(@rev_int);
         push @powerset_posns, [@rev_int];
         my @position_sets;
         foreach my $posns (@powerset_posns) {
            # due to hack in powerset code, we get a blank line instead of the full nmer
            next if (! scalar(@$posns));
            push @position_sets, [reverse(@$posns)];
         }
         foreach my $position_set (@position_sets) {
            my $nts = scalar (@$position_set);
            my $mask_a;
            my $j=0;
            for (my $i=0; $i<scalar(@$int); $i++) {
               if ( ($j < $nts) && ($$int[$i] == $$position_set[$j]) ) {
                  $$mask_a[$i] = 1;
                  $j++;
               }
               else {
                  $$mask_a[$i] = 0;
               }
            }
            if ($j != $nts) {
               print "WARNING:: mask count j does not equal number of positions ... incomplete mask", "\n";
               print $j, " != ", $nts, "\n", br;
            }
            my $seen_seq_h;
            for (my $count=0; $count<scalar (@$sequences_a); $count++) {
               my $seq_substr;
               foreach my $seq_pos (@$position_set) {
                  $seq_substr .= substr($$sequences_a[$count], ($seq_pos-1), 1);
               }
               $$seen_seq_h{$seq_substr}++;
            }
            my $name = join($JOIN_CHAR, (map {$_ + $INDEX_OFFSET} @$position_set));

            my $full_seq_pos = -1;
            if (! (exists $$WM{"name2pos"}{$full_nts}{$full_name})) {
               print "WARNING :: NAME2POS does not exist in full matrix calculation\n";
            }
            else {
               $full_seq_pos = $$WM{"name2pos"}{$full_nts}{$full_name};
            }
       
            my $subencoding_hm = return_subencoding_hmatrix($nts);
            my @subencoding_hm_keys = keys %$subencoding_hm;
            my $encode_line = scalar(@{$$subencoding_hm{$subencoding_hm_keys[0]}});
            my $seq2int_h = return_seq2int($nts);

            my $num_usable_cols = $encode_line;
            my @posset_coeff_names = @{memoized_return_coefficient_names(@$position_set)};
            my @use_thetas = splice(@$thetas,0,$num_usable_cols);
            for (my $idx=0; $idx<scalar(@use_thetas); $idx++) {
               if ($use_interaction) {
                  $use_thetas[$idx] = $$coeffName_h{$posset_coeff_names[$idx]};
               }
            }
#         print STDERR "Testing:: Usable columns $num_usable_cols out of ", $encode_line,"\n";
            foreach my $seq (sort {$$full_seq2int_h{$a} <=> $$full_seq2int_h{$b}} keys %$full_seq2int_h) {
               my $full_row = $$full_seq2int_h{$seq};
               my @seq_a = split //, $seq;
               if (scalar(@seq_a) != scalar(@$mask_a)) {
                  print "WARNING :: Sequence not size of mask", "\n";
               }
               my $sub_seq = '';
               for (my $i=0; $i<scalar(@seq_a); $i++) {
                 $sub_seq .= $seq_a[$i]    if ($$mask_a[$i] == 1);    
               }
               my $row = $$seq2int_h{$sub_seq};
               if ( (exists $$seen_seq_h{$sub_seq}) && ($use_interaction) ) {
                  my @use_encodings = @{$$subencoding_hm{$row}};
                  if (scalar(@use_encodings) != scalar(@use_thetas)) {
                     print STDERR "ERROR:: Number of dependant encoding variables (", scalar(@use_encodings), ") does not equal the number of thetas(", scalar(@use_thetas), ") that were pulled out\n";
                  }
                  else {
                     $$WM{$WM_type}{$full_nts}[$full_row][$full_seq_pos] += 0;
                     for (my $idx=0; $idx<scalar(@use_thetas); $idx++) {
                        $$WM{$WM_type}{$full_nts}[$full_row][$full_seq_pos] += $use_encodings[$idx]*$use_thetas[$idx];
                     }
                  }
               }
               else {
                  $$WM{$WM_type}{$full_nts}[$full_row][$full_seq_pos] = '-';
               }
            }
	 }
      }
#   print STDERR "DONE ...\n";
   }


   uniquify_based_on(\%$WM, "CoeffNames", "Coeffs", "Pvals");


   return $WM;
}


sub uniquify_based_on {

   my $WM = shift @_;

   my $base_type = shift @_;
   my @other_types = @_;


   my @temp_a;
   my @types_a;
   my $base_num_elem = scalar @{$$WM{$base_type}};
   push @temp_a, [@{$$WM{$base_type}}];
   push @types_a, $base_type;
   delete $$WM{$base_type};
   foreach my $type (@other_types) {
      my $num_elem = scalar @{$$WM{$type}};
      if ($num_elem != $base_num_elem) {
         print STDERR "ERROR:: Number of elements not the same in type ($type) and base_type ($base_type)\n";
      }
      else {
         push @temp_a, [@{$$WM{$type}}];
         push @types_a, $type;
         delete $$WM{$type};
      }
   }
   my %seen;
   for (my $col=0; $col < $base_num_elem; $col++) {
      next      if ($seen{$temp_a[0][$col]}++);
      for (my $row=0; $row < scalar(@temp_a); $row++) {
         push @{$$WM{$types_a[$row]}}, $temp_a[$row][$col];
      }
   }
   return;
}

#
# Finally, we have the printing functions
#


sub relativify_pwm {
   my $pwm = shift @_;
  
   my $relative_pwm;
   for (my $j=0; $j<scalar(@{$$pwm[0]}); $j++) {
      my $max = $$pwm[0][$j];
       print ($max);
      for (my $i=1; $i<scalar(@{$pwm}); $i++) {
         $max = ($$pwm[$i][$j] > $max) ? $$pwm[$i][$j] : $max;
      }
      for (my $i=0; $i<scalar(@{$pwm}); $i++) {
         $$relative_pwm[$i][$j] = $$pwm[$i][$j] - $max;
      }
   }

   return $relative_pwm;
}

sub web_print_logo {
   my $table_a = shift @_;

   my $pos2char_h = undef;
   if (@_) {
      $pos2char_h = shift @_;
   }
   my $labels;
   foreach my $pos (keys %$pos2char_h) {
      $$labels[$pos] = $$pos2char_h{$pos};
   }

   my $params;
   if ($ZOOM_LOGOS) {
      $$params{"maximize_logo_visibility"} = 1;
   }
   $$params{"svg_height"} = $LOGO_HEIGHT_CM   if (defined $LOGO_HEIGHT_CM);
   $$params{"svg_width"} = $LOGO_WIDTH_CM   if (defined $LOGO_WIDTH_CM);

   my $e_matrix = hack_energy_matrix($table_a);
   my $p_matrix = energy2probability_matrix($e_matrix, $params);
   my $svg = return_svg_logo($p_matrix, $params, $labels);
   my $base64png = svg2embeddedbase64png($svg);
   my $cmd = "convert $svg ";
   print "$svg\n";
   #print "$base64png\n";
   #print '<img alt="Embedded Image" src="data:image/png;base64,' . $base64png . '" />'. "\n"; 
}

sub web_print_elogo {
   my $table_a = shift @_;

   my $pos2char_h = undef;
   if (@_) {
      $pos2char_h = shift @_;
   }
   my $labels;
   foreach my $pos (keys %$pos2char_h) {
      $$labels[$pos] = $$pos2char_h{$pos};
   }

   my $params;
   $$params{"svg_height"} = $LOGO_HEIGHT_CM   if (defined $LOGO_HEIGHT_CM);
   $$params{"svg_width"} = $LOGO_WIDTH_CM   if (defined $LOGO_WIDTH_CM);

   my $e_matrix = hack_energy_matrix($table_a);
   my $svg = return_svg_elogo($e_matrix, $params, $labels);
   my $base64png = svg2embeddedbase64png($svg);
   print '<img alt="Embedded Image" src="data:image/png;base64,' . $base64png . '" />'. "\n"; 
}

sub hack_energy_matrix {
   my $table_a = shift @_;
   my $rows = scalar (@$table_a);
   my $cols = scalar (@{$$table_a[0]});

   my $e_matrix;
   for (my $j=0; $j<$cols; $j++) {
      my $nonzero = 0;
      my $num_nonzero = 0;
      for (my $i=0; $i<$rows; $i++) {
         next   if ($$table_a[$i][$j] eq '-');
         if ($$table_a[$i][$j] != 0) {
            $nonzero = 1;
            $num_nonzero++;
         }
      }
      for (my $i=0; $i<$rows; $i++) {
         if ($$table_a[$i][$j] eq '-') {
            $$e_matrix[$i][$j] = '-';
         }
         else {
            $$e_matrix[$i][$j] = $$table_a[$i][$j];
            if ( ($nonzero == 0) && ($num_nonzero > 1) ) {
               $$e_matrix[$i][$j] = '-';
            }
         }
      }
   }
   return $e_matrix;
}

sub svg2embeddedbase64png {
   my $svg = shift @_;
 
#   my $png = `export LD_LIBRARY_PATH=$INKSCAPE_PATH; echo \'$svg\' | $INKSCAPE_PROGRAM -z --export-png "/dev/stdout" /dev/stdin`;
   #my $png = `echo \'test\' | ./$INKSCAPE_PROGRAM --export-png "/dev/stdout" /dev/stdin`;
   open (PNG, "export LD_LIBRARY_PATH=$INKSCAPE_PATH; echo \'$svg\' | $INKSCAPE_PROGRAM -z -e \"/dev/stdout\" /dev/stdin |");
   #print "export LD_LIBRARY_PATH=$INKSCAPE_PATH; echo \'$svg\' | $INKSCAPE_PROGRAM -z -e \"/dev/stdout\" /dev/stdin |";
   binmode(PNG);
   my ($buf, $data, $n); 
   my $PNG_HEADER = '89 50 4e 47 0d ';
   my $test_string='';
   while (($n = read PNG, $data, 1) != 0) { 
      if (sprintf("%02x", ord($data)) eq '0a') {
         last    if ($test_string eq $PNG_HEADER);
         $test_string = '';
         $buf = '';
      }
      else {
         $test_string .= sprintf("%02x", ord($data)) . ' ';
         $buf .= $data;
      }
   }
   $buf .= $data;
   while (($n = read PNG, $data, 1) != 0) { 
      $buf .= $data;
   }
   close PNG;
   return encode_base64($buf);
}


sub return_transpose_matrix {
   my $matrix_a    = shift @_;
   my $pos2char_h  = shift @_;



   my $transpose_matrix = '';
   my $precision = logb(4, scalar(@{$matrix_a}));
   my @column_order = map {$_->[0]} sort {$a->[1] cmp $b->[1]} map {[$_, pack('C'.$precision, (split /$JOIN_CHAR/, $$pos2char_h{$_}))]}  keys %$pos2char_h;
   
   for (my $x=0; $x<scalar(@column_order); $x++) {
      my $j = $column_order[$x];
      $transpose_matrix .= $$pos2char_h{$column_order[$x]} . ':';     
      for (my $i=0; $i<scalar(@{$matrix_a}); $i++) {
         if ($$matrix_a[$i][$j]) {
            $transpose_matrix .= ' ' . $$matrix_a[$i][$j];
         }
         else {
            $transpose_matrix .= ' ' . '-';
         }
      }
      $transpose_matrix .= "\n";
   }
   return $transpose_matrix;
}



sub web_print_matrix {
   my $matrix_a    = shift @_;
   my $name        = shift @_;

   my $int_to_name_h = undef;
   if (@_) {
      $int_to_name_h = shift @_;
   }

   my $pos2char_h = undef;
   if (@_) {
      $pos2char_h = shift @_;
   }

   print '<table border="1">';
   print '<tr align="CENTER">';
   print td($name);

   my $precision = logb(4, scalar(@{$matrix_a}));
   my @column_order = (0 .. $#{$$matrix_a[0]});

   if (scalar (keys %$pos2char_h) == scalar(@column_order)) {
      @column_order = map {$_->[0]} sort {$a->[1] cmp $b->[1]} map {[$_, pack('C'.$precision, (split /$JOIN_CHAR/, $$pos2char_h{$_}))]}  keys %$pos2char_h;
      
      for (my $x=0; $x<scalar(@column_order); $x++) {
         print td($$pos2char_h{$column_order[$x]});
      }
   }
   else {
      for (my $x=0; $x<scalar(@column_order); $x++) {
         print td($x);
      }
   }
   print '</tr>';

   for (my $i=0; $i<scalar(@{$matrix_a}); $i++) {
      print '<tr align="CENTER">';
      
      my $base_string = '';
      if (! (exists $$int_to_name_h{$i})) {
         my $bases = dec2base($i, 4, $precision);
         foreach my $base (@$bases) {
            $base_string .= $int_to_base_h{$base};
         }
      }
      else {
         $base_string = $$int_to_name_h{$i};
      }
      print td(b($base_string));
      for (my $x=0; $x<scalar(@column_order); $x++) {
         my $j = $column_order[$x];
         my $val = round($$matrix_a[$i][$j], $PRECISION);
         print td($val);
      }
      print '</tr>';
   }
   print '</table>';
}
    
sub print_matrix {
   my $matrix_a    = shift @_;
   my $name        = shift @_;

   my $int_to_name_h = undef;
   if (@_) {
      $int_to_name_h = shift @_;
   }

   my $pos2char_h = undef;
   if (@_) {
      $pos2char_h = shift @_;
   }

   print "\n\n ------ MATRIX PRINT ------\n";
   print "#" . $name;

   my @column_order = (0 .. $#{$$matrix_a[0]});

   if (scalar (keys %$pos2char_h) == scalar(@column_order)) {
      @column_order = sort {$$pos2char_h{$a} cmp $$pos2char_h{$b}} (keys %$pos2char_h);
      for (my $x=0; $x<scalar(@column_order); $x++) {
         print "\t", $$pos2char_h{$column_order[$x]};
      }
   }
   else {
      for (my $x=0; $x<scalar(@column_order); $x++) {
         print "\t", $x;
      }
   }

   print "\n";

   my $precision = logb(4, scalar(@{$matrix_a}));
   for (my $i=0; $i<scalar(@{$matrix_a}); $i++) {
      my $base_string = '';
      if (! (exists $$int_to_name_h{$i})) {
         my $bases = dec2base($i, 4, $precision);
         foreach my $base (@$bases) {
            $base_string .= $int_to_base_h{$base};
         }
      }
      else {
         $base_string = $$int_to_name_h{$i};
         print '   ';
      }
      print $base_string . ' | ';
      for (my $x=0; $x<scalar(@column_order); $x++) {
         my $j = $column_order[$x];
         if ($$matrix_a[$i][$j] ne '-') {
            print "\t", sprintf("%1.3f", $$matrix_a[$i][$j]);
         }
         else {
            print "\t", sprintf("%5s", $$matrix_a[$i][$j]);
         }
      }
      print "\n";
   }
   print "\n";
   print "\n\n --------------------------\n";
}    

sub web_print_table {
   my $table_a = shift @_;
   
   print "\n";
   print '<table border="1">', "\n";
   foreach my $row (@$table_a) {
      print '<tr align="CENTER">', "\n";
      my @cols = @$row;
      my $name = shift @cols;
      print '<td>' . '<b>' . $name . '</b>' . '</td>';
      foreach my $col (@cols) {
         print '<td>', $col, '</td>';
      }
      print '</tr>', "\n";
   }
   print '</table>', "\n";
   print br, "\n";
}


sub print_table {
   my $table_a = shift @_;

   foreach my $row (@$table_a) {
      print join("\t", @$row), "\n";
   }
   print "\n\n";
}


sub web_print_OvE_plot {
   my $observed_a = shift @_;
   my $expected_a = shift @_;

   my $potential_ps_file = '';
   my $R_EXECUTABLE = 'R --vanilla -q --slave';
   local (*READ, *WRITE);
   my $pid = open2(\*READ, \*WRITE, $R_EXECUTABLE);
                
   # read input file
   print WRITE 'data<-read.table(stdin(),sep=",")', "\n";

   if (scalar (@$observed_a) != scalar(@$expected_a)) {
      print "WARNING :: Observed != Expected\n";
   }

   my $min = min(@$observed_a, @$expected_a);
   my $max = max(@$observed_a, @$expected_a);
   my $lim = 'c(' . $min . ',' . $max . ')';
   for (my $i=0; $i<scalar(@$observed_a); $i++) {
      print WRITE $$observed_a[$i], ',', $$expected_a[$i], "\n";
   }
   print WRITE "\n";

   my $HEIGHT = 5;
   my $WIDTH  = 5;

   print WRITE 'postscript(file=' . '""' . ', command="cat", horizontal=FALSE, paper="special", height=' . $HEIGHT .', width=' . $WIDTH . ')' . "\n";
   
   print WRITE 'attach(data)' . "\n";
   print WRITE 'plot(data$V1, data$V2, ' . ' xlim='. $lim . ', ylim=' . $lim . ',' . ' main = list("Observed vs. Expected", col="black"), xlab="Observed", ylab ="Expected", type="p", col="black")' . "\n";
   print WRITE 'abline(0,1)' . "\n";
   
   # closes current device
   print WRITE 'dev.off()' . "\n";

   # closes R process
   print WRITE 'q()' . "\n";
   while (my $ps_line = <READ>) {
      chomp $ps_line;
      $potential_ps_file .= $ps_line . "\n";
      last     if  ($ps_line eq '%%EOF');
   }
   while (my $leftover_trash = <READ>) {
      # get rid of remining junk
   }
   
   my $w = $WIDTH*72;
   my $h = $HEIGHT*72;
   open (PNG, "echo \"$potential_ps_file\" | gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -dGraphicsAlphaBits=4 -dDEVICEWIDTHPOINTS=$w -dDEVICEHEIGHTPOINTS=$h -sOutputFile=%stdout - |");
   binmode(PNG);
   my ($buf, $data, $n); 
   my $PNG_HEADER = '89 50 4e 47 0d ';
   my $test_string='';
   while (($n = read PNG, $data, 1) != 0) { 
      if (sprintf("%02x", ord($data)) eq '0a') {
        last    if ($test_string eq $PNG_HEADER);
         $test_string = '';
         $buf = '';
      }
      else {
         $test_string .= sprintf("%02x", ord($data)) . ' ';
         $buf .= $data;
      }
   }
   $buf .= $data;
   while (($n = read PNG, $data, 1) != 0) { 
      $buf .= $data;
   }
   close PNG;
   `echo \"$potential_ps_file\" > trash_potential_ps.ps`;

   my $png_webgraph_base64 = encode_base64($buf);
   close READ;
   close WRITE;
   waitpid($pid, 0);
   print '<img alt="Embedded R Graphs" src="data:image/png;base64,' . $png_webgraph_base64 . '" />';
}


sub max {
   my $max = $_[0];

   foreach (@_)  { $max = ($max > $_) ? $max : $_; }

   return $max;
}

sub min {
   my $min = $_[0];

   foreach (@_)  { $min = ($min < $_) ? $min : $_; }

   return $min;
}


