/*
 * Copyright (c) 2012 LabKey Corporation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.github.discvrseq.util;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CigarUtil;

import java.util.Iterator;

/**
 * User: bbimber
 * Date: 8/9/12
 * Time: 9:34 PM
 */

/**
 * Iterator that transverses the elements of a CIGAR string, accumulating information at each position of the alignment.
 *
 */
public class CigarPositionIterable implements Iterable<CigarPositionIterable.PositionInfo>
{
    public static final byte INDEL_CHARACTER = (byte)'-';
    public static final byte AMBIGUITY_CHARACTER = (byte)'N';

    private SAMRecord _record;

    public CigarPositionIterable(SAMRecord record)
    {
        _record = record;
    }

    public CigarIterator iterator()
    {
        return new CigarIterator(this);
    }

    public class CigarIterator implements Iterator<PositionInfo>
    {
        private SAMRecord _record;

        private Integer[] _readPositions;
        private Integer[] _refPositions;
        private char[] _explodedCigar;
        private int _pos = 0;

        /**
         * Prepare to iterate the CIGAR string of this record
         * @param iterable A CigarPositionIterable instance
         */
        public CigarIterator(CigarPositionIterable iterable)
        {
            _record = iterable._record;
            initializeCigar();
        }

        private void initializeCigar()
        {
            Cigar c = _record.getCigar();
            _explodedCigar = CigarUtil.cigarArrayFromString(c.toString());

            int readPos = 0; //0-based
            int refPos = _record.getAlignmentStart() - 1; //0-based

            _readPositions = new Integer[_explodedCigar.length];
            _refPositions = new Integer[_explodedCigar.length];
            int i = 0;
            for (char el : _explodedCigar)
            {
                CigarOperator op = CigarOperator.valueOf(Character.toString(el));
                if (op.consumesReadBases())
                {
                    _readPositions[i] = readPos;
                    readPos++;
                }
                else
                {
                    _readPositions[i] = -1;
                }

                if (op.consumesReferenceBases())
                {
                    _refPositions[i] = refPos;
                    refPos++;
                }
                else
                {
                    _refPositions[i] = -1;
                }
                i++;
            }
        }

        public void remove()
        {
            throw new UnsupportedOperationException("Can not remove elements in CIGAR via an iterator!");
        }

        /**
         * @return true if there are more positions of the CIGAR string to iterate
         */
        public boolean hasNext()
        {
            return _pos < _readPositions.length;
        }

        /**
         * @return The next PositionInfo in the iteration
         */
        public PositionInfo next()
        {
            if (_pos >= _readPositions.length)
                return null;

            PositionInfo info = new PositionInfo(_record, _pos, _explodedCigar, _readPositions, _refPositions);
            _pos++;
            return info;
        }
    }

    /**
     * Describes a specific position in an alignment, including the position relative to the start of both the reference and read
     * sequences.
     */
    public class PositionInfo
    {
        private SAMRecord _record;
        private CigarOperator _op;
        private int _pos;
        private int _readPos;
        private int _refPos;
        private int _indel = 0;

        private int _lastReadPos;
        private int _lastRefPos;

        public PositionInfo(SAMRecord record, int pos, char[] ops, Integer[] readPos, Integer[] refPos)
        {
            _record = record;
            _pos = pos;
            _op = CigarOperator.valueOf(Character.toString(ops[pos]));
            _readPos = readPos[pos];
            _refPos = refPos[pos];

            if (_readPos > -1)
                _lastReadPos = _readPos;
            else
            {
                int i = _pos;
                while (i >= 0)
                {
                    if (readPos[i] > -1)
                    {
                        _lastReadPos = readPos[i];
                        _indel = i - _pos;
                        break;
                    }
                    i--;
                }
            }

            if (_refPos > -1)
                _lastRefPos = _refPos;
            else
            {
                int i = _pos;
                while (i >= 0)
                {
                    if (refPos[i] > -1)
                    {
                        _lastRefPos = refPos[i];
                        _indel = _pos - i;
                        break;
                    }
                    i--;
                }
            }
        }

        /**
         * @return The zero-based position relative to the start of the reference, -1 indicates an insertion
         */
        public int getRefPosition()
        {
            return _refPos;
        }

        /**
         * @return The zero-based position relative to the start of the read, -1 indicates a deletion
         */
        public int getReadPosition()
        {
            return _readPos;
        }

        /**
         * @return This will return the length of the indel at this position: 0 for no indel, positive for an insertion relative to the reference,
         * negative for a deletion relative to the reference.  For example, the 1st base of an insertion is 1, the second is 2, etc.
         * The first base of a deletion is -1, second is -2.  This is patterned after the Bio-SamTools (perl) Bio::DB::Bam::Pileup::indel method.
         */
        public int getIndel()
        {
            return _indel;
        }

        /**
         * @return The length of the insertion at this position.  For non-indels or deletions, it will be zero.  The first position of an insertion is 0 (ie. the base overlapping with the references), the second is 1, etc.
         */
        public int getInsertIndex()
        {
            return getIndel() <= 0 ? 0 : getIndel();
        }

        /**
         * @return The last read position that overlapped the reference, using zero-based coordinates.  This is primarily used for deletions relative to
         * the reference, since those return -1 as the read position.  For non-indels, this will return the same value as getReadPosition()
         */
        public int getLastReadPosition()
        {
            return _lastReadPos;
        }

        /**
         * @return Returns the last reference position that overlapped the read, using zero-based coordinates.  Primarily used for insertions relative
         * to the reference, in order to find the previous reference position.  For non-indels, this will return the same value as getRefPosition()
         */
        public int getLastRefPosition()
        {
            return _lastRefPos;
        }

        /**
         * @param referenceBases An array representing the sequence of the reference
         * @return The reference base at this position.  '-' indicates an insertion.
         */
        public byte getReferenceBase(byte[] referenceBases)
        {
            return isInsertion() ? INDEL_CHARACTER : referenceBases[_refPos];
        }

        /**
         * @return The read base at this position.  '-' indicates a deletion.
         */
        public byte getReadBase()
        {
            return isDel() ? INDEL_CHARACTER : _record.getReadBases()[_readPos];
        }

        /**
         * @return True if this position is skipped, meaning the CIGAR operator is a soft clip (S) or a skipped region (N)
         */
        public boolean isSkipped()
        {
            return _op.equals(CigarOperator.SOFT_CLIP) || _op.equals(CigarOperator.HARD_CLIP) || _op.equals(CigarOperator.SKIPPED_REGION);
        }

        /**
         * @return True if this position is an indel (insertion or deletion) relative to the reference
         */
        public boolean isIndel()
        {
            return isDel() || isInsertion();
        }

        /**
         * @return True if this position is an insertion relative to the reference
         */
        public boolean isInsertion()
        {
            return _op.equals(CigarOperator.INSERTION);
        }

        /**
         * @return True if this position is a deletion relative to the reference
         */
        public boolean isDel()
        {
            return _op.equals(CigarOperator.DELETION);
        }

        /**
         * @return True if this position overlaps the reference, meaning it consumes both read and references bases
         */
        public boolean overlapsReference()
        {
            return _op.consumesReadBases() && _op.consumesReferenceBases();
        }

        public boolean includeInSnpCount()
        {
            return isIndel() || (overlapsReference() && !isSkipped());
        }

        public String getReferenceName()
        {
            return _record.getReferenceName();
        }

        /**
         * @return The CigarOperation at this position
         */
        public CigarOperator getCigarOperator()
        {
            return _op;
        }

        /**
         * @return The SAMRecord associated with this alignment position
         */
        public SAMRecord getRecord()
        {
            return _record;
        }

        /**
         * @return The base quality at this position
         */
        public int getBaseQuality()
        {
            return _record.getBaseQualities()[getLastReadPosition()];
        }
    }
}

