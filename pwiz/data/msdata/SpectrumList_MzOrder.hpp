//
// $Id$
//
//
// Original author: Brian Pratt <bspratt .at. proteinms.net>
// AI assistance: Claude Code (Claude Opus 5) <noreply .at. anthropic.com>
//
// Copyright 2026 University of Washington - Seattle, WA
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//


#ifndef _SPECTRUMLIST_MZORDER_HPP_
#define _SPECTRUMLIST_MZORDER_HPP_


#include "pwiz/utility/misc/Export.hpp"
#include "SpectrumListWrapper.hpp"
#include <atomic>


namespace pwiz {
namespace msdata {


/// Presents peaks in ascending m/z order whatever order the writer used.
///
/// Ascending m/z is nowhere required by the mzML specification, but it is what every consumer
/// assumes: extraction binary searches the m/z axis, so a spectrum presented in any other order
/// makes the search land nowhere useful and the chromatogram comes out empty, with no error at
/// all. Writers that present some other order do exist - one shipped peaks in ascending
/// intensity - and msconvert passes that order through untouched, so a file converted from one
/// is broken for everything downstream. This is applied on read so that no consumer of the
/// library has to know about it.
///
/// The question is settled from the first few spectra of a file rather than re-asked for every
/// spectrum, since walking every m/z array of every file to catch a rare writer is a cost the
/// whole world would pay for the few. The first spectrum alone will not do: early scans can
/// precede the sample and carry almost no peaks, and a spectrum with two of them ascends half
/// the time by chance. So the two verdicts are not symmetric - one spectrum out of order proves
/// the writer does not sort however few peaks it holds, while peaks found in order are only
/// believed from a spectrum with enough of them to mean it, and until one arrives the checking
/// continues.
class PWIZ_API_DECL SpectrumList_MzOrder : public SpectrumListWrapper
{
    public:

    SpectrumList_MzOrder(const SpectrumListPtr& inner);

    virtual SpectrumPtr spectrum(size_t index, bool getBinaryData = false) const;
    virtual SpectrumPtr spectrum(size_t index, DetailLevel detailLevel) const;
    virtual bool benefitsFromWorkerThreads() const;

    /// True once a spectrum has actually been put back in order, which also means the file was
    /// not written in ascending m/z order.
    bool reorderedAnySpectrum() const {return reordered_.load();}

    private:

    /// What the file has shown so far about the way its writer orders peaks.
    enum class Verdict {unsettled, writerSortsByMz, writerDoesNotSortByMz};

    void ensureMzAscending(const SpectrumPtr& spectrum) const;

    // Mutable and atomic because spectrum() is const and is called from worker threads.
    mutable std::atomic<Verdict> verdict_;
    mutable std::atomic<bool> reordered_;
};


} // namespace msdata
} // namespace pwiz


#endif // _SPECTRUMLIST_MZORDER_HPP_
