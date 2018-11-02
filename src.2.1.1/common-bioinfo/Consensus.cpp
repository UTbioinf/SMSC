#include "Consensus.h"

namespace loon
{// namespace loon starts

/*static*/ Aligner Consensus::aligner;

Consensus::Consensus():
    left_delimiter(INF), right_delimiter(INF)
{
}

void Consensus::swap(Consensus& t)
{
    std::swap(left_delimiter, t.left_delimiter);
    std::swap(right_delimiter, t.right_delimiter);
    std::vector<Nucleotide>::swap(t);
}

void Consensus::swap(std::vector<Nucleotide>& t)
{
    std::vector<Nucleotide>::swap(t);
}

void Consensus::set_str(const std::string& str)
{
    size_t n = str.length();
    this->resize( n );
    for(size_t i = 0; i < n; ++i)
        this->at(i).init_base(str[i], true);
    left_delimiter = INF;
    right_delimiter = INF;
}

std::string& Consensus::toString()
{
    str.clear();
    str.reserve(this->size());
    for(std::vector<Nucleotide>::iterator it = this->begin(); it != this->end(); ++it)
        str.push_back( it->get_max_base() );
    return str;
}

void Consensus::toString(std::string& str)
{
    str.clear();
    str.reserve(this->size());
    for(std::vector<Nucleotide>::iterator it = this->begin(); it != this->end(); ++it)
        str.push_back( it->get_max_base() );
}

void Consensus::consensus_init(const std::string& qry)
{
    aligner.set_strs( toString(), qry );
    aligner.set_cmp_obj();
}

void Consensus::consensus_error(const std::string& qry)
{
    ERROR_PRINT("consensus failed, lol, save data!!!");
    save_binary("err_ref.consensus.bin");
    toString();
    save_fasta("err_ref", str, "err_ref.fasta");
    save_fasta("err_qry", qry, "err_qry.fasta");
    exit(1);
}

void Consensus::consensus_result(const std::string& qry)
{
    DEBUG_PRINT("Consensus::Consensus_result()");
    DEBUG_PRINT("len: (XX, YY) = (%lu, %lu)", aligner.XX->length(), aligner.YY->length());
    DEBUG_PRINT("XX = (%lu, %lu), YY = (%lu, %lu), indels.size() = %lu",
            aligner.XX_start, aligner.XX_end,
            aligner.YY_start, aligner.YY_end,
            aligner.indels.size());

    size_t pos_x = aligner.XX_start - 1;
    size_t pos_y = aligner.YY_start - 1;
    bool left_delimiter_set = false;
    bool right_delimiter_set = false;

    // compute aligned part
    std::vector<Nucleotide> middle_part;
    middle_part.reserve( aligner.XX_end - aligner.XX_start + 1 );

    //long discard_cnt = 0, insert_cnt = 0;

    for(std::vector<long>::const_iterator indel_it = aligner.indels.begin();
            indel_it != aligner.indels.end(); ++indel_it)
    {
        long steps = abs(*indel_it);
        while(--steps > 0)
        {
            this->at( pos_x ).set_base( qry[ pos_y ] );
            middle_part.push_back( this->at( pos_x ) );
            ++pos_x;
            ++pos_y;
        }

        if(*indel_it > 0)
        {
            this->at( pos_x ).set_base( '-' );
            if(!this->at( pos_x ).is_dash())
                middle_part.push_back( this->at(pos_x) );
            //else
            //    ++discard_cnt;
            ++pos_x;
        }
        else
        {
            //++insert_cnt;
            middle_part.push_back( Nucleotide(qry[pos_y]) );
            ++pos_y;
        }
    }

    //DEBUG_PRINT("discard_cnt = %ld, insert_cnt = %ld", discard_cnt, insert_cnt);

    assert(aligner.XX_end - pos_x == aligner.YY_end - pos_y);

    while(pos_x < aligner.XX_end)
    {
        this->at( pos_x ).set_base( qry[pos_y] );
        middle_part.push_back( this->at( pos_x ) );
        ++pos_x;
        ++pos_y;
    }

    assert( pos_x == aligner.XX_end );
    assert( pos_y == aligner.YY_end );

    // compute left padding
    std::vector<Nucleotide> left_padding;
    if(aligner.XX_start == 1)
    {
        /*
        ......##########
        ################
        */
        for(size_t i = 1; i < aligner.YY_start; ++i)
            left_padding.push_back( Nucleotide(qry[i-1], true) );

        // left_delimiter
        // right_delimiter
        right_delimiter = 0;
        right_delimiter_set = true;
    }
    else if(aligner.YY_start == 1)
    {
        /*
        ##############
        ......########
        */
        left_padding.assign(this->begin(), this->begin() + aligner.XX_start);
    }
    else
    {
        ERROR_PRINT("left is not totally aligned");
        exit(1);
    }

    // compute right padding
    std::vector<Nucleotide> right_padding;
    pos_x = aligner.XX_end;
    pos_y = aligner.YY_end;
    if(aligner.XX_end == aligner.XX->length())
    {
        /*
        #########............
        #####################
        */
        while(pos_y < aligner.YY->length())
            right_padding.push_back( Nucleotide( qry[pos_y++], true ) );
        // left_delimiter
        left_delimiter = left_padding.size() + middle_part.size() + right_padding.size() - 1;
        left_delimiter_set = true;
    }
    else if(aligner.YY_end == aligner.YY->length())
    {
        /*
        #######################
        ########...............
        */
        right_padding.assign(this->begin() + pos_x, this->end());
    }
    else
    {
        ERROR_PRINT("right is not totally aligned");
        exit(1);
    }

    DEBUG_PRINT("left_padding = %lu, middle_part = %lu, right_padding = %lu", left_padding.size(), middle_part.size(), right_padding.size());

    // merge
    this->clear();
    this->swap( left_padding );

    if(!right_delimiter_set)
        right_delimiter = this->size();

    this->insert( this->end(), middle_part.begin(), middle_part.end() );
    
    if(!left_delimiter_set)
        left_delimiter = this->size() - 1;

    this->insert( this->end(), right_padding.begin(), right_padding.end() );
}

void Consensus::consensus_full_strict(const std::string& str, bool unaligned_stop/* = false */)
{
    consensus_init(str);
    //aligner.set_parameters(20, 5);
    if(aligner.align_full( false ))
        consensus_result(str);
    else if(unaligned_stop)
        consensus_error(str);
    else
        throw aligner.error_code;
}

void Consensus::consensus_full_strict(void (*compute_deltas_func)(const std::string&,
                const std::string&, nucmer_Deltas&, size_t, int),
        const std::string& str, bool unaligned_stop/* = false */)
{
    consensus_init(str);
    //aligner.set_parameters(20, 5);
    if(aligner.align_full( compute_deltas_func, false ))
        consensus_result(str);
    else if(unaligned_stop)
        consensus_error(str);
    else
        throw aligner.error_code;
}

void Consensus::consensus_full(const std::string& str, bool unaligned_stop/* = false */)
{
    consensus_init(str);
    if(aligner.align_full( true ))
        consensus_result(str);
    else if(unaligned_stop)
        consensus_error(str);
    else
        throw aligner.error_code;
}

void Consensus::consensus_full(void (*compute_deltas_func)(const std::string&,
                const std::string&, nucmer_Deltas&, size_t, int),
        const std::string& str, bool unaligned_stop/* = false */)
{
    consensus_init(str);
    if(aligner.align_full( compute_deltas_func, true ))
        consensus_result(str);
    else if(unaligned_stop)
        consensus_error(str);
    else
        throw aligner.error_code;
}


void Consensus::consensus_left(const std::string& str, size_t right_most/* = INF */,
        bool unaligned_stop/* = false*/)
{
    consensus_init( str );
    if(aligner.align_left(true, right_most))
        consensus_result( str );
    else if(unaligned_stop)
        consensus_error( str );
    else
        throw aligner.error_code;
}

void Consensus::consensus_left(void (*compute_deltas_func)(const std::string&,
                const std::string&, nucmer_Deltas&, size_t, int),
        const std::string& str,
        size_t right_most/* = INF */, bool unaligned_stop/* = false*/)
{
    consensus_init( str );
    if(aligner.align_left(compute_deltas_func, true, right_most))
        consensus_result( str );
    else if(unaligned_stop)
        consensus_error( str );
    else
        throw aligner.error_code;
}

void Consensus::consensus_right(const std::string& str, size_t left_most/* = INF */,
        bool unaligned_stop/* = false */)
{
    consensus_init( str );
    //aligner.set_parameters(20, 5);
    if(aligner.align_right(true, left_most))
        consensus_result( str );
    else if(unaligned_stop)
        consensus_error( str );
    else
        throw aligner.error_code;
}

void Consensus::consensus_right(void (*compute_deltas_func)(const std::string&,
                const std::string&, nucmer_Deltas&, size_t, int),
        const std::string& str,
        size_t left_most/* = INF */, bool unaligned_stop/* = false */)
{
    consensus_init( str );
    //aligner.set_parameters(20, 5);
    if(aligner.align_right(compute_deltas_func, true, left_most))
        consensus_result( str );
    else if(unaligned_stop)
        consensus_error( str );
    else
        throw aligner.error_code;
}

void Consensus::consensus_left_auto_delimiter(const std::string& str,
        bool unaligned_stop/* = false */)
{
    if(left_delimiter == INF)
        left_delimiter = this->size() - 1;
    consensus_left(str, left_delimiter, unaligned_stop);
}

void Consensus::consensus_left_auto_delimiter(void (*compute_deltas_func)(const std::string&,
                const std::string&, nucmer_Deltas&, size_t, int),
        const std::string& str, bool unaligned_stop/* = false */)
{
    if(left_delimiter == INF)
        left_delimiter = this->size() - 1;
    consensus_left(compute_deltas_func, str, left_delimiter, unaligned_stop);
}

void Consensus::consensus_right_auto_delimiter(const std::string& str,
        bool unaligned_stop/* = false */)
{
#ifdef DEBUG
    size_t tmp = right_delimiter;
    size_t tmp_ref_length = this->size();
    size_t tmp_qry_length = str.length();
#endif

    if(right_delimiter == INF)
        right_delimiter = 0;
    consensus_right(str, right_delimiter, unaligned_stop);

    DEBUG_PRINT("last delimiter: %lu, read length = %lu, total length = %lu, current delimiter: %lu", 
            tmp, str.length(), this->size(), right_delimiter);
}

void Consensus::consensus_right_auto_delimiter(void (*compute_deltas_func)(const std::string&,
                const std::string&, nucmer_Deltas&, size_t, int),
        const std::string& str, bool unaligned_stop/* = false */)
{
#ifdef DEBUG
    size_t tmp = right_delimiter;
    size_t tmp_ref_length = this->size();
    size_t tmp_qry_length = str.length();
#endif

    if(right_delimiter == INF)
        right_delimiter = 0;
    consensus_right(compute_deltas_func, str, right_delimiter, unaligned_stop);

    DEBUG_PRINT("last delimiter: %lu, read length = %lu, total length = %lu, current delimiter: %lu", 
            tmp, str.length(), this->size(), right_delimiter);
}

void Consensus::save_binary(std::ostream& out)
{

    size_t n = this->size();
    out.write(reinterpret_cast<const char*>(&n), sizeof(size_t));
    for(size_t i = 0; i < n; ++i)
    {
        out.write(reinterpret_cast<const char*>(this->at(i).bases_cnt), sizeof(this->at(i).bases_cnt));
        out.write(reinterpret_cast<const char*>(&(this->at(i).max_cnt)), sizeof(this->at(i).max_cnt));
        out.write(reinterpret_cast<const char*>(&(this->at(i).total_cnt)), sizeof(this->at(i).total_cnt));
        out.write(&(this->at(i).ch), sizeof(char));
    }
}

void Consensus::save_binary(const std::string& fname)
{
    std::ofstream fout(fname.c_str(), std::ofstream::binary);
    check_file_open(fout, fname);

    save_binary(fout);

    fout.close();
}

void Consensus::read_binary(std::istream& in)
{
    size_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof(size_t));
    this->reserve(n);

    for(size_t i = 0; i < n; ++i)
    {
        this->push_back(Nucleotide());
        Nucleotide& tmp = this->back();

        in.read(reinterpret_cast<char*>(tmp.bases_cnt), sizeof(tmp.bases_cnt));
        in.read(reinterpret_cast<char*>(&tmp.max_cnt), sizeof(tmp.max_cnt));
        in.read(reinterpret_cast<char*>(&tmp.total_cnt), sizeof(tmp.total_cnt));
        in.read(&tmp.ch, sizeof(char));
    }
}

void Consensus::read_binary(const std::string& fname)
{
    std::ifstream fin(fname.c_str(), std::ifstream::binary);
    check_file_open(fin, fname);

    read_binary(fin);

    fin.close();
}

}// namespace loon ends
